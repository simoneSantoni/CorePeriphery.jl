#!/usr/bin/env bash
set -euo pipefail

output_tsv="${1:-benchmark/versioned_performance.tsv}"
output_md="${2:-benchmark/VERSIONED_PERFORMANCE.md}"
mode="${3:---full}"
if [[ "${mode}" != "--full" && "${mode}" != "--quick" ]]; then
    printf 'third argument must be --full or --quick\n' >&2
    exit 2
fi
temporary_directory="$(mktemp -d)"
trap 'rm -rf "${temporary_directory}"' EXIT

source_revision="$(git rev-parse HEAD)"
if [[ -n "$(git status --porcelain)" ]]; then
    source_revision="${source_revision}+working-tree"
fi
printf 'julia_version\tsource_revision\tthreads\tblas_threads\tcpu\trepetitions\tbenchmark\tn\tstorage\tmedian_ms\tminimum_ms\tmaximum_ms\tbytes\n' > "${output_tsv}"
for channel in 1.10 1.12; do
    raw="${temporary_directory}/julia-${channel}.tsv"
    benchmark_arguments=(--check)
    if [[ "${mode}" == "--quick" ]]; then
        benchmark_arguments=(--quick --check)
    fi
    julia +"${channel}" --threads=2 --project=. --startup-file=no \
        benchmark/scaling.jl "${benchmark_arguments[@]}" > "${raw}"
    context="$(awk -F '\t' '$1 == "context" {print}' "${raw}")"
    threads="$(printf '%s\n' "${context}" | awk -F '\t' '{for(i=1;i<=NF;i++) if($i ~ /^threads=/){sub(/^threads=/,"",$i); print $i}}')"
    blas_threads="$(printf '%s\n' "${context}" | awk -F '\t' '{for(i=1;i<=NF;i++) if($i ~ /^blas_threads=/){sub(/^blas_threads=/,"",$i); print $i}}')"
    cpu="$(printf '%s\n' "${context}" | awk -F '\t' '{for(i=1;i<=NF;i++) if($i ~ /^cpu=/){sub(/^cpu=/,"",$i); print $i}}')"
    repetitions="$(printf '%s\n' "${context}" | awk -F '\t' '{for(i=1;i<=NF;i++) if($i ~ /^repetitions=/){sub(/^repetitions=/,"",$i); print $i}}')"
    awk -F '\t' -v version="${channel}" -v source_revision="${source_revision}" \
        -v threads="${threads}" -v blas_threads="${blas_threads}" \
        -v cpu="${cpu}" -v repetitions="${repetitions}" \
        'BEGIN {OFS="\t"} $1 != "benchmark" && NF == 7 {
            print version,source_revision,threads,blas_threads,cpu,repetitions,$1,$2,$3,$4,$5,$6,$7
        }' \
        "${raw}" >> "${output_tsv}"
done

python3 benchmark/summarize_versioned_performance.py "${output_tsv}" "${output_md}"
printf 'wrote %s and %s\n' "${output_tsv}" "${output_md}"
