#!/usr/bin/env python3
# tcga_kirc_prep.py
# Upstream data cleaning for TCGA-KIRC. Converts raw TCGA downloads into
# the standard expression matrix and clinical table shared by all
# downstream figures and tables.
#
# Input directory (default: <repo>/data/raw; override with --data-dir or
# TCGA_KIRC_RAW):
#   TCGA-KIRC_gene_annotation.csv
#   mrna/TCGA-KIRC_tpm_mrna.csv
#   TCGA-KIRC_clinical_survival.csv
#   TCGA-KIRC_clinical_info.csv
#
# Output directory (default: <repo>/data/tcga_kirc; override with
# --output-dir or NCC_DATA_ROOT):
#   KIRC_expr_log2_tpm.csv       log2(TPM+1) matrix
#   KIRC_expr_tpm.csv            raw TPM matrix
#   KIRC_clinical_cleaned.csv    cleaned clinical table
#   KIRC_data_summary.json       filtering summary
#
# Four cleaning steps:
#   1. keep protein_coding genes only
#   2. keep primary tumour samples only (sample_code == '01')
#   3. keep the earliest vial per patient (min vial letter, 01A before 01B)
#   4. low-expression filter: mean TPM >= 1.0 AND expressed in >= 20% samples
# Step 4 is an extra step and is easy to omit from Methods descriptions.

import sys, io, os, csv, re, argparse, json, math
from collections import Counter, defaultdict

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                                  errors='replace')

MIN_MEAN_TPM = 1.0
MIN_SAMPLES_EXPRESSED = 0.2

BARCODE_RE = re.compile(r'^TCGA-\w{2}-\w{4}-(\d{2})([A-C])(?:-.*)?$')


def _find_project_root(start):
    cur = start
    for _ in range(5):
        if os.path.isdir(os.path.join(cur, "data")):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    return start


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT = _find_project_root(SCRIPT_DIR)

DATA_DIR = os.path.join(PROJ_ROOT, "data", "raw")
OUTPUT_DIR = os.path.join(PROJ_ROOT, "data", "tcga_kirc")


def parse_args():
    p = argparse.ArgumentParser(
        description="TCGA-KIRC preprocessing: protein-coding + primary tumour "
                    "+ earliest vial + low-expression filter")
    p.add_argument("--data-dir",
                   default=os.environ.get("TCGA_KIRC_RAW", DATA_DIR),
                   help="raw TCGA-KIRC directory")
    default_out = os.environ.get("NCC_DATA_ROOT")
    default_out = (os.path.join(default_out, "tcga_kirc")
                   if default_out else OUTPUT_DIR)
    p.add_argument("--output-dir", default=default_out,
                   help="output directory for cleaned data")
    p.add_argument("--min-mean-tpm", type=float, default=MIN_MEAN_TPM)
    p.add_argument("--min-samples-expressed", type=float,
                   default=MIN_SAMPLES_EXPRESSED)
    return p.parse_args()


def parse_barcode(barcode):
    parts = barcode.split('-')
    if len(parts) >= 4:
        patient = '-'.join(parts[:3])
        sample_code = parts[3]
        sample_type_code = sample_code[:2]
        vial = sample_code[2:] if len(sample_code) > 2 else ''
        type_map = {
            '01': 'Primary Solid Tumor',
            '02': 'Recurrent Solid Tumor',
            '03': 'Primary Blood Derived Cancer',
            '05': 'Additional New Primary',
            '06': 'Metastatic',
            '11': 'Solid Tissue Normal',
        }
        sample_type = type_map.get(sample_type_code,
                                   f'Unknown({sample_type_code})')
        return patient, sample_type_code, sample_type, vial
    return barcode, '??', 'Unknown', ''


def load_gene_annotation():
    print("\n[1/6] Loading gene annotation...")
    gene_types = {}
    with open(f"{DATA_DIR}/TCGA-KIRC_gene_annotation.csv",
              'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gene_types[row['gene_name']] = row['gene_type']
    protein_coding = {g for g, t in gene_types.items()
                      if t == 'protein_coding'}
    print(f"      total genes: {len(gene_types)}")
    print(f"      protein_coding: {len(protein_coding)}")
    return protein_coding, gene_types


def load_and_clean_expression(protein_coding_genes):
    print("\n[2/6] Loading mRNA TPM matrix...")
    expr_file = f"{DATA_DIR}/mrna/TCGA-KIRC_tpm_mrna.csv"
    gene_names = []

    with open(expr_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        samples_raw = header[1:]

        sample_info = []
        for s in samples_raw:
            patient, code, stype, vial = parse_barcode(s)
            sample_info.append({
                'full_barcode': s,
                'patient': patient,
                'sample_code': code,
                'sample_type': stype,
                'vial': vial,
            })

        gene_expr = []
        for row in reader:
            gene = row[0]
            if gene not in protein_coding_genes:
                continue
            try:
                expr = [float(x) if x else 0.0 for x in row[1:]]
            except (ValueError, IndexError):
                continue
            gene_names.append(gene)
            gene_expr.append(expr)

    print(f"      loaded: {len(gene_names)} protein_coding genes")
    print(f"      samples: {len(samples_raw)}")

    type_counts = Counter(s['sample_type'] for s in sample_info)
    print("\n      sample type distribution:")
    for t, c in type_counts.most_common():
        print(f"        {t}: {c}")

    tumor_indices = [i for i, s in enumerate(sample_info)
                     if s['sample_code'] == '01']
    tumor_samples = [sample_info[i] for i in tumor_indices]
    tumor_expr = [[gene_expr[g][i] for i in tumor_indices]
                  for g in range(len(gene_names))]
    print(f"\n      primary tumour samples (01A/01B): {len(tumor_samples)}")

    patient_samples = defaultdict(list)
    for i, s in enumerate(tumor_samples):
        patient_samples[s['patient']].append(i)

    dup_patients = {p: idxs for p, idxs in patient_samples.items()
                    if len(idxs) > 1}
    print(f"      patients with multiple tumour samples: {len(dup_patients)}")

    keep_indices = []
    dropped = 0
    for p, idxs in sorted(patient_samples.items()):
        best = min(idxs, key=lambda i: tumor_samples[i]['vial'])
        keep_indices.append(best)
        if len(idxs) > 1:
            dropped += len(idxs) - 1

    tumor_samples_final = [tumor_samples[i] for i in keep_indices]
    tumor_expr_final = [[tumor_expr[g][i] for i in keep_indices]
                        for g in range(len(gene_names))]
    print(f"      after deduplication: {len(tumor_samples_final)} samples "
          f"({dropped} duplicates dropped)")

    print(f"\n      low-expression filter (mean TPM >= {MIN_MEAN_TPM}, "
          f"expressed in >= {MIN_SAMPLES_EXPRESSED*100:.0f}% samples)...")
    n_samples = len(tumor_samples_final)
    keep_genes = []
    keep_expr = []
    for gene, expr_row in zip(gene_names, tumor_expr_final):
        mean_tpm = sum(expr_row) / n_samples
        frac_expressed = sum(1 for v in expr_row if v > 0) / n_samples
        if mean_tpm >= MIN_MEAN_TPM and frac_expressed >= MIN_SAMPLES_EXPRESSED:
            keep_genes.append(gene)
            keep_expr.append(expr_row)
    print(f"      genes before: {len(gene_names)} -> after: {len(keep_genes)}")

    print("\n      applying log2(TPM+1)...")
    log_expr = [[math.log2(v + 1) for v in row] for row in keep_expr]

    sample_ids = [s['patient'] for s in tumor_samples_final]

    return {
        'genes': keep_genes,
        'samples': sample_ids,
        'full_barcodes': [s['full_barcode'] for s in tumor_samples_final],
        'expr_tpm': keep_expr,
        'expr_log2': log_expr,
        'n_genes': len(keep_genes),
        'n_samples': len(sample_ids),
        'all_sample_info': sample_info,
        'tumor_sample_info': tumor_samples_final,
    }


def load_and_clean_clinical(tumor_patients):
    print("\n[3/6] Loading clinical data...")

    survival = {}
    with open(f"{DATA_DIR}/TCGA-KIRC_clinical_survival.csv",
              'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            pid = row['case_submitter_id']
            survival[pid] = {
                'OS': row.get('OS', ''),
                'OS.time': row.get('OS.time', ''),
                'gender': row.get('gender', ''),
                'age': row.get('age', ''),
                'T': row.get('T', ''),
                'M': row.get('M', ''),
                'N': row.get('N', ''),
                'stage': row.get('stage', ''),
            }
    print(f"      survival records: {len(survival)}")

    clinical_detail = {}
    with open(f"{DATA_DIR}/TCGA-KIRC_clinical_info.csv",
              'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            pid = row.get('submitter_id', '')
            if pid:
                clinical_detail[pid] = {
                    'ajcc_pathologic_stage': row.get('ajcc_pathologic_stage', ''),
                    'ajcc_pathologic_t': row.get('ajcc_pathologic_t', ''),
                    'ajcc_pathologic_n': row.get('ajcc_pathologic_n', ''),
                    'ajcc_pathologic_m': row.get('ajcc_pathologic_m', ''),
                    'tumor_grade': row.get('tumor_grade', ''),
                    'laterality': row.get('laterality', ''),
                    'age_at_diagnosis': row.get('age_at_diagnosis', ''),
                    'gender': row.get('gender', ''),
                    'race': row.get('race', ''),
                    'ethnicity': row.get('ethnicity', ''),
                    'morphology': row.get('morphology', ''),
                    'tissue_or_organ_of_origin': row.get('tissue_or_organ_of_origin', ''),
                    'primary_diagnosis': row.get('primary_diagnosis', ''),
                    'prior_malignancy': row.get('prior_malignancy', ''),
                    'prior_treatment': row.get('prior_treatment', ''),
                    'vital_status': row.get('vital_status', ''),
                    'days_to_death': row.get('days_to_death', ''),
                    'days_to_last_follow_up': row.get('days_to_last_follow_up', ''),
                    'year_of_diagnosis': row.get('year_of_diagnosis', ''),
                    'site_of_resection_or_biopsy': row.get('site_of_resection_or_biopsy', ''),
                    'tumor_focality': row.get('tumor_focality', ''),
                }
    print(f"      detailed clinical records: {len(clinical_detail)}")

    print("\n[4/6] Merging and cleaning clinical phenotypes...")
    merged = []
    missing_surv = 0
    missing_detail = 0

    for pid in sorted(tumor_patients):
        surv = survival.get(pid, {})
        detail = clinical_detail.get(pid, {})
        if not surv:
            missing_surv += 1
        if not detail:
            missing_detail += 1

        age_raw = surv.get('age', '') or detail.get('age_at_diagnosis', '')
        try:
            age_val = float(age_raw)
            if age_val > 1000:
                age_val = round(age_val / 365.25)
            elif age_val > 100:
                age_val = None
        except (ValueError, TypeError):
            age_val = None

        os_raw = surv.get('OS', '')
        os_time_raw = surv.get('OS.time', '')
        try:
            os = int(os_raw) if os_raw not in ('', 'NA', None) else None
        except (ValueError, TypeError):
            os = None
        try:
            os_time = float(os_time_raw) if os_time_raw not in ('', 'NA', None) else None
        except (ValueError, TypeError):
            os_time = None

        if detail.get('days_to_death', ''):
            try:
                dtd = float(detail['days_to_death'])
                if dtd > 0:
                    os = 1
                    os_time = dtd
            except (ValueError, TypeError):
                pass
        elif detail.get('days_to_last_follow_up', '') and os_time is None:
            try:
                dlf = float(detail['days_to_last_follow_up'])
                if dlf > 0:
                    os = 0
                    os_time = dlf
            except (ValueError, TypeError):
                pass

        stage_raw = surv.get('stage', '') or detail.get('ajcc_pathologic_stage', '')
        stage_clean = stage_raw.replace('Stage ', '').strip() if stage_raw else ''

        grade_raw = detail.get('tumor_grade', '')
        grade_clean = ''
        if grade_raw:
            m = re.search(r'G(\d)', grade_raw)
            if m:
                grade_clean = m.group(1)

        t_raw = surv.get('T', '') or detail.get('ajcc_pathologic_t', '')
        n_raw = surv.get('N', '') or detail.get('ajcc_pathologic_n', '')
        m_raw = surv.get('M', '') or detail.get('ajcc_pathologic_m', '')

        def clean_tnm(val):
            if not val or val in ('', 'NA', 'NX', 'Not Reported'):
                return ''
            val = val.strip()
            for prefix in ['T', 'N', 'M']:
                if val.startswith(prefix) and len(val) > 1 and val[1].isdigit():
                    val = val[1:]
            return val

        gender_raw = surv.get('gender', '') or detail.get('gender', '')
        gender_clean = gender_raw.strip().upper() if gender_raw else ''
        if gender_clean in ('MALE', 'M'):
            gender_clean = 'MALE'
        elif gender_clean in ('FEMALE', 'F'):
            gender_clean = 'FEMALE'

        laterality_raw = detail.get('laterality', '')
        laterality_clean = laterality_raw.strip().title() if laterality_raw else ''

        race_raw = detail.get('race', '')
        race_clean = race_raw.strip() if race_raw else ''

        merged.append({
            'patient_id': pid,
            'age': age_val,
            'gender': gender_clean,
            'race': race_clean,
            'laterality': laterality_clean,
            'os': os,
            'os_time_days': os_time,
            'os_time_years': round(os_time / 365.25, 2) if os_time else None,
            'stage': stage_clean,
            'grade': grade_clean,
            'T': clean_tnm(t_raw),
            'N': clean_tnm(n_raw),
            'M': clean_tnm(m_raw),
            'vital_status': detail.get('vital_status', ''),
            'year_of_diagnosis': detail.get('year_of_diagnosis', ''),
            'prior_malignancy': detail.get('prior_malignancy', ''),
        })

    if missing_surv:
        print(f"      WARNING missing survival data: {missing_surv}")
    if missing_detail:
        print(f"      WARNING missing detailed clinical data: {missing_detail}")

    n_total = len(merged)
    n_os_avail = sum(1 for m in merged if m['os'] is not None)
    n_death = sum(1 for m in merged if m['os'] == 1)
    n_alive = sum(1 for m in merged if m['os'] == 0)
    n_age_avail = sum(1 for m in merged if m['age'] is not None)
    n_stage_avail = sum(1 for m in merged if m['stage'])
    n_grade_avail = sum(1 for m in merged if m['grade'])

    print("\n      cleaning complete:")
    print(f"        total patients: {n_total}")
    print(f"        OS available: {n_os_avail} (dead={n_death}, alive={n_alive})")
    print(f"        age available: {n_age_avail}")
    print(f"        stage available: {n_stage_avail}")
    print(f"        grade available: {n_grade_avail}")

    stage_dist = Counter(m['stage'] for m in merged if m['stage'])
    print("\n      stage distribution:")
    for s in ['I', 'II', 'III', 'IV']:
        print(f"        Stage {s}: {stage_dist.get(s, 0)}")

    grade_dist = Counter(m['grade'] for m in merged if m['grade'])
    print("\n      grade distribution:")
    for g in ['1', '2', '3', '4']:
        print(f"        G{g}: {grade_dist.get(g, 0)}")

    return merged


def match_and_merge(expr_data, clinical_data):
    print("\n[5/6] Matching expression matrix with clinical data...")

    expr_patients = set(expr_data['samples'])
    clin_patients = set(c['patient_id'] for c in clinical_data)
    common = expr_patients & clin_patients
    expr_only = expr_patients - clin_patients
    clin_only = clin_patients - expr_patients

    print(f"      expression patients: {len(expr_patients)}")
    print(f"      clinical patients:   {len(clin_patients)}")
    print(f"      overlap:             {len(common)}")
    if expr_only:
        print(f"      WARNING expression-only: {len(expr_only)}")
    if clin_only:
        print(f"      WARNING clinical-only:   {len(clin_only)}")

    patient_to_expr_idx = {}
    for i, pid in enumerate(expr_data['samples']):
        if pid in common:
            patient_to_expr_idx[pid] = i
    patient_to_clin = {c['patient_id']: c for c in clinical_data
                       if c['patient_id'] in common}

    common_sorted = sorted(common)
    expr_indices = [patient_to_expr_idx[pid] for pid in common_sorted]
    clin_records = [patient_to_clin[pid] for pid in common_sorted]

    final_genes = expr_data['genes']
    final_expr_log2 = [[expr_data['expr_log2'][g][i] for i in expr_indices]
                       for g in range(len(final_genes))]
    final_expr_tpm = [[expr_data['expr_tpm'][g][i] for i in expr_indices]
                      for g in range(len(final_genes))]

    print(f"\n      final dataset: {len(common_sorted)} patients × "
          f"{len(final_genes)} genes")

    return {
        'patients': common_sorted,
        'genes': final_genes,
        'expr_log2': final_expr_log2,
        'expr_tpm': final_expr_tpm,
        'clinical': clin_records,
        'n_patients': len(common_sorted),
        'n_genes': len(final_genes),
    }


def save_cleaned_data(final_data):
    print("\n[6/6] Saving cleaned data...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    expr_file = f"{OUTPUT_DIR}/KIRC_expr_log2_tpm.csv"
    print(f"      expression matrix (log2 TPM+1): {expr_file}")
    with open(expr_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['gene'] + final_data['patients'])
        for g_idx, gene in enumerate(final_data['genes']):
            writer.writerow([gene] + final_data['expr_log2'][g_idx])

    tpm_file = f"{OUTPUT_DIR}/KIRC_expr_tpm.csv"
    print(f"      expression matrix (TPM):        {tpm_file}")
    with open(tpm_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['gene'] + final_data['patients'])
        for g_idx, gene in enumerate(final_data['genes']):
            writer.writerow([gene] + final_data['expr_tpm'][g_idx])

    clin_file = f"{OUTPUT_DIR}/KIRC_clinical_cleaned.csv"
    print(f"      clinical table:                 {clin_file}")
    with open(clin_file, 'w', newline='', encoding='utf-8') as f:
        fields = list(final_data['clinical'][0].keys())
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for rec in final_data['clinical']:
            writer.writerow({k: '' if v is None else v for k, v in rec.items()})

    summary = {
        'dataset': 'TCGA-KIRC',
        'n_patients': final_data['n_patients'],
        'n_genes': final_data['n_genes'],
        'gene_type': 'protein_coding',
        'expression_type': 'log2(TPM+1)',
        'min_mean_tpm': MIN_MEAN_TPM,
        'min_samples_expressed': MIN_SAMPLES_EXPRESSED,
        'cleaning_steps': [
            'keep protein_coding genes',
            'keep primary tumor samples (sample_code == 01)',
            'keep earliest vial per patient (min vial letter, 01A before 01B)',
            'low-expression filter',
        ],
        'clinical_fields': list(final_data['clinical'][0].keys()) if final_data['clinical'] else [],
        'os_events': sum(1 for c in final_data['clinical'] if c['os'] == 1),
        'os_censored': sum(1 for c in final_data['clinical'] if c['os'] == 0),
        'stage_distribution': dict(Counter(c['stage'] for c in final_data['clinical'] if c['stage'])),
        'grade_distribution': dict(Counter(c['grade'] for c in final_data['clinical'] if c['grade'])),
    }
    summary_file = f"{OUTPUT_DIR}/KIRC_data_summary.json"
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(f"      data summary:                   {summary_file}")

    print("=" * 60)
    print("                    cleaning summary")
    print("=" * 60)
    print(f"  dataset:   TCGA-KIRC")
    print(f"  matrix:    {final_data['n_genes']} genes × "
          f"{final_data['n_patients']} samples")
    print(f"  log2 TPM:  {expr_file}")
    print(f"  clinical:  {clin_file}")
    print(f"  OS events: {summary['os_events']} (dead) / "
          f"{summary['os_censored']} (censored)")
    print("=" * 60)


def main():
    global DATA_DIR, OUTPUT_DIR, MIN_MEAN_TPM, MIN_SAMPLES_EXPRESSED

    args = parse_args()
    DATA_DIR = args.data_dir
    OUTPUT_DIR = args.output_dir
    MIN_MEAN_TPM = args.min_mean_tpm
    MIN_SAMPLES_EXPRESSED = args.min_samples_expressed

    print("=" * 60)
    print("  TCGA-KIRC data cleaning pipeline")
    print("=" * 60)
    print(f"  DATA_DIR   = {DATA_DIR}")
    print(f"  OUTPUT_DIR = {OUTPUT_DIR}")
    print(f"  low-expression filter: mean TPM >= {MIN_MEAN_TPM}, "
          f"expressed in >= {MIN_SAMPLES_EXPRESSED*100:.0f}% samples")

    protein_coding, _ = load_gene_annotation()
    expr_data = load_and_clean_expression(protein_coding)
    tumor_patients = set(expr_data['samples'])
    clinical_data = load_and_clean_clinical(tumor_patients)
    final_data = match_and_merge(expr_data, clinical_data)
    save_cleaned_data(final_data)

    print("\n  Done.")


if __name__ == '__main__':
    main()