---
license: mit
private: true
configs:
- config_name: default
  data_files:
  - split: train
    path: data/train-*
dataset_info:
  features:
  - name: model_id
    dtype: string
  - name: model_name
    dtype: string
  - name: provider
    dtype: string
  - name: model_type
    dtype: string
  - name: parameters_billions
    dtype: float64
  - name: license
    dtype: string
  - name: context_window
    dtype: int64
  - name: modality
    dtype: string
  - name: architecture
    dtype: string
  - name: aime2026_score
    dtype: float64
  - name: evasionBench_score
    dtype: float64
  - name: gpqa_score
    dtype: float64
  - name: gsm8k_score
    dtype: float64
  - name: hle_score
    dtype: float64
  - name: hmmt2026_score
    dtype: float64
  - name: mmluPro_score
    dtype: float64
  - name: olmOcr_score
    dtype: float64
  - name: swePro_score
    dtype: float64
  - name: sweVerified_score
    dtype: float64
  - name: terminalBench_score
    dtype: float64
  - name: aggregate_score
    dtype: float64
  - name: coverage_count
    dtype: int64
  - name: coverage_percent
    dtype: float64
  splits:
  - name: train
    num_bytes: 25680
    num_examples: 105
  download_size: 17894
  dataset_size: 25680
---

