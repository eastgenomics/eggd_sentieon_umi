# eggd_sentieon_umi: Technical Walkthrough

*2026-02-20T16:55:32Z by Showboat 0.6.0*
<!-- showboat-id: 5083b1d9-9a38-43e2-aee5-c96c0b1cb9a7 -->

## Overview

`eggd_sentieon_umi` is a DNAnexus applet (v1.0.0) that implements a UMI-aware short-read alignment pipeline using Sentieon's high-performance genomics toolkit. It takes raw FASTQ files as input and produces a sorted, duplicate-removed BAM or CRAM file, leveraging Unique Molecular Identifiers (UMIs) to distinguish PCR duplicates from true biological reads.

The applet runs on DNAnexus under Ubuntu 20.04, targets a 36-core instance (`mem1_ssd1_v2_x36` in AWS eu-central-1), and ships the Sentieon suite (v202503.01) and samtools as bundled resources.

## Repository Structure

<details>
<summary>Repository file listing</summary>

```bash
find . -not -path '*/.git/*' -not -name 'technical_walkthrough.md' -type f | sed 's|^\./||' | sort
```

```output
dxapp.json
Readme.md
resources/home/dnanexus/helper_funcs.sh
resources/home/dnanexus/license_auth.sh
resources/home/dnanexus/license_setup.sh
resources/home/dnanexus/README.txt
resources/usr/bin/samtools
resources/usr/local/sentieon-genomics-202503.01/bin/sentieon
resources/usr/local/sentieon-genomics-202503.01/doc/LEGAL
resources/usr/local/sentieon-genomics-202503.01/doc/licsrvr.service
resources/usr/local/sentieon-genomics-202503.01/doc/licsrvr.sh
resources/usr/local/sentieon-genomics-202503.01/doc/RELNOTES
resources/usr/local/sentieon-genomics-202503.01/libexec/bamslice
resources/usr/local/sentieon-genomics-202503.01/libexec/bwa
resources/usr/local/sentieon-genomics-202503.01/libexec/driver
resources/usr/local/sentieon-genomics-202503.01/libexec/fqidx
resources/usr/local/sentieon-genomics-202503.01/libexec/licclnt
resources/usr/local/sentieon-genomics-202503.01/libexec/licsrvr
resources/usr/local/sentieon-genomics-202503.01/libexec/LongReadUtil
resources/usr/local/sentieon-genomics-202503.01/libexec/minimap2
resources/usr/local/sentieon-genomics-202503.01/libexec/plot
resources/usr/local/sentieon-genomics-202503.01/libexec/pyexec
resources/usr/local/sentieon-genomics-202503.01/libexec/rcat
resources/usr/local/sentieon-genomics-202503.01/libexec/STAR
resources/usr/local/sentieon-genomics-202503.01/libexec/umi
resources/usr/local/sentieon-genomics-202503.01/libexec/util
resources/usr/local/sentieon-genomics-202503.01/lib/libblas.so.3.2.1
resources/usr/local/sentieon-genomics-202503.01/lib/libbz2.so.1.0.4
resources/usr/local/sentieon-genomics-202503.01/lib/libgfortran.so.3.0.0
resources/usr/local/sentieon-genomics-202503.01/lib/libgomp.so.1.0.0
resources/usr/local/sentieon-genomics-202503.01/lib/libhdf5.so.6.0.4
resources/usr/local/sentieon-genomics-202503.01/lib/liblapack.so.3.2.1
resources/usr/local/sentieon-genomics-202503.01/lib/liblzma.so.0.0.0
resources/usr/local/sentieon-genomics-202503.01/lib/libz.so.1.2.8
resources/usr/local/sentieon-genomics-202503.01/lib/python/sentieon/vcflib/bgzf.py
resources/usr/local/sentieon-genomics-202503.01/lib/python/sentieon/vcflib/compat.py
resources/usr/local/sentieon-genomics-202503.01/lib/python/sentieon/vcflib/__init__.py
resources/usr/local/sentieon-genomics-202503.01/lib/python/sentieon/vcflib/sharder.py
resources/usr/local/sentieon-genomics-202503.01/lib/python/sentieon/vcflib/tabix.py
resources/usr/local/sentieon-genomics-202503.01/lib/python/sentieon/vcflib/tribble.py
resources/usr/local/sentieon-genomics-202503.01/lib/python/sentieon/vcflib/vcf.py
resources/usr/local/sentieon-genomics-202503.01/share/funcs
src/code.sh
```

</details>

The top-level layout is:

- **`dxapp.json`** — DNAnexus app manifest: input/output spec, runtime configuration, instance requirements, and access rules.
- **`src/code.sh`** — The sole bash entry-point executed by the DNAnexus worker.
- **`resources/home/dnanexus/`** — Helper scripts placed into `/home/dnanexus/` on the worker: licence discovery, authentication, and refresh logic.
- **`resources/usr/bin/samtools`** — Bundled samtools binary.
- **`resources/usr/local/sentieon-genomics-202503.01/`** — The full Sentieon v202503.01 distribution: the `sentieon` CLI, libexec workers (`bwa`, `umi`, `driver`, `bamslice`, `minimap2`, `STAR`, …), shared libraries, and Python VCF utilities.

## App Manifest — `dxapp.json`

The manifest defines the DNAnexus app contract. Key sections are shown below.

<details>
<summary>App identity</summary>

```bash
jq '{name,title,summary,version,categories,openSource}' dxapp.json
```

```output
{
  "name": "eggd_sentieon_umi",
  "title": "eggd_sentieon_umi",
  "summary": "Performs pre-procesing of UMI FASTQ files through UMI tag extraction and barcode-aware duplicate removal and consensus calling, and final alignment to create a sorted and deduplicated BAM file from FASTQ files",
  "version": "1.0.0",
  "categories": [
    "Read Mapping",
    "DNAseq"
  ],
  "openSource": false
}
```

</details>

### Input specification

The app defines 17 inputs across four logical groups.

<details>
<summary>Full input spec</summary>

```bash
jq -r '.inputSpec[] | [.name, .class, (.optional // false | tostring), (.group // "(top-level)"), .label] | @tsv' dxapp.json | column -t -s $'\t'
```

```output
readsidx_fastqgzs     array:file  true   (top-level)         Reads (UMI sequence)
reads_fastqgzs        array:file  false  (top-level)         Reads
reads2_fastqgzs       array:file  false  (top-level)         Reads (right mates)
rg_info_csv           file        true   (top-level)         Read group information
genomebwaindex_targz  file        false  (top-level)         BWA reference genome index
genome_fastagz        file        false  (top-level)         Reference genome FASTA file
read_group_platform   string      false  Sample Information  Read group platform
sample                string      true   Sample Information  Sample ID
stream_inputs         boolean     false  Input Options       Stream Input Reads?
output_format         string      false  Output Options      Output a BAM or a CRAM file?
output_metrics        boolean     false  Output Options      Calculate and output the BAM metrics?
output_md5sum         boolean     false  Output Options      Calculate output md5sum
read_template         string      false  Algorithm Options   UMI read template
duplex_umi            boolean     false  Algorithm Options   Duplex UMI
extra_bwa_options     string      true   Algorithm Options   Extra BWA Options
bwa_nonverbose        boolean     false  Advanced Options    Reduce BWA logs
bam_compression       string      false  Advanced Options    Intermediate BAM compression
```

</details>

Inputs fall into four groups:

| Group | Inputs | Purpose |
|---|---|---|
| *(top-level)* | `reads_fastqgzs`, `reads2_fastqgzs`, `readsidx_fastqgzs`, `rg_info_csv`, `genomebwaindex_targz`, `genome_fastagz` | Core data: FASTQ files, optional UMI index reads, optional read-group CSV, reference genome and BWA index |
| Sample Information | `read_group_platform`, `sample` | BAM `@RG` tags; sample name defaults to FASTQ filename prefix |
| Input / Output Options | `stream_inputs`, `output_format`, `output_metrics`, `output_md5sum` | Stream FASTQs via `dx cat`; choose BAM vs CRAM; enable metrics PDFs; emit md5 checksums |
| Algorithm / Advanced | `read_template`, `duplex_umi`, `extra_bwa_options`, `bwa_nonverbose`, `bam_compression` | UMI read structure syntax, duplex mode, extra BWA flags, log verbosity, intermediate compression level |

The `read_template` field uses a `<number><operator>` syntax (e.g. `8M12S+T,+T`) where `M` = molecular barcode, `S` = spacer to skip, `T` = template sequence, `+` = rest of read.

### Output specification

<details>
<summary>Full output spec</summary>

```bash
jq -r '.outputSpec[] | [.name, .class, (.optional // false | tostring), .label] | @tsv' dxapp.json | column -t -s $'\t'
```

```output
mappings_bam      file        false  Sorted mappings
mappings_bam_bai  file        false  Sorted mappings index
metrics           array:file  true   Reads stats metrics
md5sums           array:file  true   MD5 sum files for most outputs
```

</details>

### Runtime configuration

<details>
<summary>Runtime config JSON</summary>

```bash
jq '{runSpec: {interpreter,file,distribution,release,version,timeoutPolicy,execDepends,assetDepends}, regionalOptions, access, authorizedUsers}' dxapp.json
```

```output
{
  "runSpec": {
    "interpreter": null,
    "file": null,
    "distribution": null,
    "release": null,
    "version": "1.0.0",
    "timeoutPolicy": null,
    "execDepends": null,
    "assetDepends": null
  },
  "regionalOptions": {
    "aws:eu-central-1": {
      "systemRequirements": {
        "main": {
          "instanceType": "mem1_ssd1_v2_x36"
        }
      }
    }
  },
  "access": {
    "network": [
      "*"
    ],
    "allProjects": "CONTRIBUTE"
  },
  "authorizedUsers": [
    "org-emee_1"
  ]
}
```

```bash
jq '{runSpec: .runSpec, regionalOptions, access, authorizedUsers}' dxapp.json
```

```output
{
  "runSpec": {
    "timeoutPolicy": {
      "*": {
        "hours": 24
      }
    },
    "assetDepends": [
      {
        "name": "htslib_suite_asset",
        "project": "project-Fkb6Gkj433GVVvj73J7x8KbV",
        "folder": "/app_assets/htslib/htslib_v1.22.0",
        "version": "1.22"
      }
    ],
    "execDepends": [
      {
        "name": "rename"
      }
    ],
    "distribution": "Ubuntu",
    "release": "20.04",
    "version": "0",
    "file": "src/code.sh",
    "interpreter": "bash"
  },
  "regionalOptions": {
    "aws:eu-central-1": {
      "systemRequirements": {
        "main": {
          "instanceType": "mem1_ssd1_v2_x36"
        }
      }
    }
  },
  "access": {
    "network": [
      "*"
    ],
    "allProjects": "CONTRIBUTE"
  },
  "authorizedUsers": [
    "org-emee_1"
  ]
}
```

</details>

Key runtime details:

- **OS**: Ubuntu 20.04
- **Entry point**: `src/code.sh` (bash interpreter)
- **Instance**: `mem1_ssd1_v2_x36` (36 vCPUs, 72 GB RAM) in `aws:eu-central-1`
- **Timeout**: 24 hours
- **System package dependency**: `rename` (Perl rename utility, installed at job start via `apt-get`)
- **DNAnexus asset dependency**: `htslib_suite_asset` v1.22 (provides htslib/bgzf shared libraries at `/opt/dnanexus/assets/htslib/lib`)
- **Network access**: unrestricted (required for Sentieon licence server contact)
- **Project access**: `CONTRIBUTE` on all projects (required for reading licence records)

## Pipeline Walkthrough — `src/code.sh`

The entry-point script is structured as a set of modular bash functions, called sequentially at the bottom of the file. The two main modules are `module_download_inputs` and `module_umi_align`.

### Step 0 — Initialisation

Before any module runs, `code.sh` does three things:

1. Sets `set -eo pipefail` — the script exits immediately on any non-zero exit code or failed pipe.
2. Installs `gnuplot` (used later by `sentieon plot metrics` to render the PDF metrics report).
3. Sources `license_setup.sh`, which discovers the Sentieon licence, exports all required `SENTIEON_*` environment variables, and validates the licence with a `licclnt ping` before any compute is attempted.

<details>
<summary>Initialisation block — <code>src/code.sh</code> lines 1–10</summary>

```bash
sed -n '1,10p' src/code.sh
```

```output
#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -eo pipefail
apt-get update && apt-get install -y gnuplot
# Set up options
source /home/dnanexus/license_setup.sh

export LD_LIBRARY_PATH=/opt/dnanexus/assets/htslib/lib:$LD_LIBRARY_PATH
```

</details>

### Step 1 — `module_download_inputs`

This module downloads all input files in parallel, with two optimisations:

- **Streaming the genome FASTA**: The reference FASTA is never fully written to disk before indexing begins. Instead, a named FIFO (`genome_temp.fa`) lets `samtools faidx` index the stream concurrently while the full file is written to `genome/`.
- **Streaming FASTQ inputs** (when `stream_inputs=true`): FASTQ files are excluded from the bulk download and instead streamed on-the-fly with `dx cat` directly into the alignment pipe, saving the time and storage of a full download.

After unpacking, the module validates that the BWA index prefix matches the FASTA filename, and infers the sample name from the FASTQ filename if not explicitly provided.

<details>
<summary><code>module_download_inputs</code> — <code>src/code.sh</code> lines 125–187</summary>

```bash
sed -n '125,187p' src/code.sh
```

```output
module_download_inputs()
{
	# Download all inputs 
	EXCEPT_LIST="$EXCEPT_LIST --except genome_fastagz --except mappings_bam"
	dx-download-all-inputs --parallel $EXCEPT_LIST

    input_command="mv -f --"
    if [[ "$in_files" != "" ]]; then
        eval $input_command "${in_files_path[@]}" .
    fi

	# Stream and unpack FASTA genome while using samtools to index it
	wait_pids_arr=()
	mkdir genome
	genome_file="genome/${genome_fastagz_name%.gz}"
	mkfifo genome_temp.fa
	dx download "${genome_fastagz}" -o - | gunzip -c | tee genome_temp.fa > "${genome_file}" &
	gunzipGenomeFastagz_pid=$!
	wait_pids_arr+=(${gunzipGenomeFastagz_pid})
	samtools faidx genome_temp.fa &
	samtools_genome_faidx_pid=$!
	wait_pids_arr+=(${samtools_genome_faidx_pid})
	if [ "$genomebwaindex_targz" != "" ]; then
		tar --no-same-owner -zxvf "$genomebwaindex_targz_path" -C genome &
		untarGenomeIndexTargz_pid=$!
		wait_pids_arr+=(${untarGenomeIndexTargz_pid})
	fi
	
	# Wait until untar commands are done
	wait_pids "${wait_pids_arr[@]}"
	unset wait_pids_arr
	wait_pids_arr=()
	
	# Rename genome file
	genome_file_prefix="${genomebwaindex_targz_name%.bwa-index.tar.gz}"
	rename "s/genome\.fa/${genome_file_prefix}\.fa/" ./genome/genome* || true
	mv genome_temp.fa.fai "${genome_file}.fai"
	rm genome_temp.fa
	
	#check if the BWA indices and the FASTA match
	if [ "$genomebwaindex_targz" != "" ]; then
		found_right=false
		for bwt_file in genome/${genome_file_prefix}.fa*.bwt; do
			if [ "${bwt_file%.bwt}" == "${genome_file}" ]; then
				found_right=true
			fi
		done
		if [ "$found_right" == "false" ]; then
			error_report "The BWA reference genome index does not correspond to the Reference genome FASTA file. The app will exit now."
		fi
	fi
	
	# Set sample name if not given
	if [ "$sample" = "" ]; then
		sample="${reads_fastqgzs_prefix[0]}"
		sample="$(echo $sample|sed 's|_00[0-9]$||'|sed 's|_R1$||'|sed 's|_1$||')"
	else
		sample="${sample// /_}"
	fi
	if [ "$sample" = "" ]; then
		sample="sample"
	fi
}
```

</details>

### Step 2 — `module_umi_align`

This is the core alignment module. It loops over each FASTQ pair and runs a streaming pipeline for each:

<details>
<summary>Alignment pipe — <code>src/code.sh</code> lines 277–310</summary>

```bash
sed -n '277,310p' src/code.sh
```

```output
		if [ "$stream_inputs" = "true" ]; then
			{ (${SENTIEON_INSTALL_DIR}/bin/sentieon umi extract $duplex_arg ${read_template} $read_idx_path \
				<(stream_input "${reads_fastqgzs[$i]}") <(stream_input "${reads2_fastqgzs[$i]}") \
				| ${SENTIEON_INSTALL_DIR}/bin/sentieon bwa mem -t ${sentieon_procs} $extra_bwa_options \
				"${genome_file}" ${bwa_options} \
				-R "@RG\tID:${read_group_id}\tPL:${read_group_platform}\tPU:${read_group_pu}\tLB:${read_group_library}\tSM:${CASE_SAMPLE}" \
				-p -C - || echo -n 'error' ) 2>&3 \
				| ${SENTIEON_INSTALL_DIR}/bin/sentieon util sort -i - --sam2bam -o sorted_${i}.bam \
				-t ${sentieon_procs} ${util_sort_options} --bam_compression $bam_compression --block_size 2G ;}\
				3>&1 1>&2 | grep -v --line-buffered "^\[M::mem_pestat" \
				| grep -v --line-buffered "^\[M::process"| ([ "$bwa_nonverbose" != "true" ] \
				&& cat || grep -v --line-buffered "^\[M::mem_process" )
		else
			{ (${SENTIEON_INSTALL_DIR}/bin/sentieon umi extract $duplex_arg ${read_template} $read_idx_path \
				"${reads_fastqgzs_path[$i]}" "${reads2_fastqgzs_path[$i]}" \
				| ${SENTIEON_INSTALL_DIR}/bin/sentieon bwa mem -t ${sentieon_procs} $extra_bwa_options \
				"${genome_file}" ${bwa_options} \
				-R "@RG\\tID:${read_group_id}\\tPL:${read_group_platform}\\tPU:${read_group_pu}\\tLB:${read_group_library}\\tSM:${CASE_SAMPLE}" \
				-p -C - || echo -n 'error' ) 2>&3 \
				| ${SENTIEON_INSTALL_DIR}/bin/sentieon util sort -i - --sam2bam -o sorted_${i}.bam \
				-t ${sentieon_procs} ${util_sort_options} --bam_compression $bam_compression --block_size 2G ;}\
				3>&1 1>&2 | grep -v --line-buffered "^\[M::mem_pestat" \
				| grep -v --line-buffered "^\[M::process"| ([ "$bwa_nonverbose" != "true" ] \
				&& cat || grep -v --line-buffered "^\[M::mem_process" )
			rm ${reads_fastqgzs_path[$i]}
			if [ ${#reads2_fastqgzs_path[@]} -ne 0 ]; then
				rm ${reads2_fastqgzs_path[$i]}
			fi
		fi
		sorted_bam_args+=("-i")
		sorted_bam_args+=("sorted_${i}.bam")
		number_fastq=$((number_fastq - 1))
		echo "$number_fastq remaining FASTQ files to be mapped"
	done
```

</details>

The per-pair pipe (streaming mode) is:

```output
sentieon umi extract [--duplex] <read_template> R1.fq.gz R2.fq.gz
  | sentieon bwa mem -t <nproc> -p -C -R <@RG> <genome>
  | sentieon util sort --sam2bam -o sorted_N.bam
```

Three stages run end-to-end without intermediate files:

1. **`sentieon umi extract`** — Parses the read template string to locate the UMI bases in each read, strips them from the sequence, and writes them into SAM auxiliary tags (`XR` for the UMI, and others). The `-d` / `--duplex` flag enables duplex-UMI mode where both strands carry an identical UMI structure.
2. **`sentieon bwa mem`** — Multi-threaded BWA-MEM alignment. The `-p` flag indicates interleaved paired-end input; `-C` passes through SAM auxiliary tags (including the UMI tags from step 1) to the output. The `@RG` header is injected here (ID, PL, PU, LB, SM).
3. **`sentieon util sort`** — Converts SAM → BAM and sorts by coordinate, using configurable compression (`bam_compression` 1/6/9) and 2 GB block size for I/O efficiency.

File descriptor tricks (`2>&3 3>&1 1>&2`) swap stdout and stderr so that BWA progress messages (on stderr) can be selectively filtered by grep without breaking the stdout SAM stream.

### Step 3 — LocusCollector consensus + Dedup

<details>
<summary>LocusCollector + Dedup commands — <code>src/code.sh</code> lines 335–344</summary>

```bash
sed -n '335,344p' src/code.sh
```

```output

	${SENTIEON_INSTALL_DIR}/bin/sentieon driver -t ${sentieon_procs} -r "${genome_file}" \
		"${sorted_bam_args[@]}" --algo LocusCollector --consensus --umi_tag XR \
		--fun score_info score.txt.gz  "${metrics_algo_args[@]}"
	${SENTIEON_INSTALL_DIR}/bin/sentieon driver -t ${sentieon_procs} -r "${genome_file}" \
		"${sorted_bam_args[@]}" --algo Dedup --score_info score.txt.gz --metrics \
		"${CASE_PREFIX}dedup_metrics.txt" "${dedup_bam_name}.${dedup_bam_extension}"

	rm sorted_*.bam || true
}
```

</details>

After all per-pair BAMs are sorted, two passes over all BAMs run back-to-back:

**Pass 1 — `LocusCollector` (score collection)**

```output
sentieon driver -t <nproc> -r <genome> -i sorted_0.bam [-i sorted_1.bam ...]
  --algo LocusCollector --consensus --umi_tag XR --fun score_info score.txt.gz
  [--algo CoverageMetrics ... --algo GCBias ... ...]   # only if output_metrics=true
```

`LocusCollector` with `--consensus` and `--umi_tag XR` groups reads by genomic position and UMI (`XR` tag written by `umi extract`), scores each read group, and writes a compressed score file (`score.txt.gz`). Optional quality metrics (GC bias, insert size, alignment stats, quality distribution, coverage) are computed in the same single pass when `output_metrics=true`.

**Pass 2 — `Dedup` (consensus BAM)**

```output
sentieon driver -t <nproc> -r <genome> -i sorted_*.bam
  --algo Dedup --score_info score.txt.gz
  --metrics <sample>.dedup_metrics.txt
  <output>.bam   # or .cram
```

`Dedup` reads the score file and the sorted BAMs, collapses each UMI family into a single consensus read, and writes the final deduplicated output in BAM or CRAM format (controlled by `output_format`). Deduplication metrics (duplication rate, estimated library size, etc.) are always written.

Intermediate per-pair sorted BAMs are deleted after this step to free disk space.

### Step 4 — Output upload

Outputs are uploaded asynchronously using a custom `upload()` / `streaming_upload()` pattern that:

1. Moves each result file into the DNAnexus output staging directory (`~/out/<handle>/`).
2. Calls `dx upload` in the background and records the resulting file ID.
3. Optionally computes and uploads an `.md5` checksum file alongside the main output (when `output_md5sum=true`).
4. Waits for all background uploads to complete via `wait_uploads()` before registering outputs with `dx-jobutil-add-output`.
5. Calls `dx-upload-all-outputs --parallel` for any remaining files not individually tracked (primarily the metrics directory).

The `wait_uploads()` function polls for sentinel files (`uploaded_<handle>.file`) created by `dx upload` on success, and surfaces upload errors via `uploaded_<handle>.file.error` files without immediately aborting the job.

## Licence Management — `resources/home/dnanexus/`

Sentieon requires a valid licence to run. The three helper scripts handle discovery, setup, and ongoing refresh of the licence token.

<details>
<summary><code>license_setup.sh</code></summary>

```bash
cat resources/home/dnanexus/license_setup.sh
```

```output
#!/bin/bash

source ~/helper_funcs.sh
SENTIEON_INSTALL_DIR=$(sentieon_install_dir)
sentieon_procs=$(nproc)
# Setup license
# find corresponding project where license token information is and download it
find_license_project
if [ "$license_project" == "" ]; then
	error_report "You do not have a license to use this app. Please contact Sentieon at info@sentieon.com to request a license. The app will exit now."
fi
download_license_token

# Kick token refresh script
bash ~/license_auth.sh $license_project&

# Set license settings
job_tag=$(cat ~/Sentieon_License|jq '.job_tag')
license_server_location=$(cat ~/Sentieon_License|jq '.license_server_location'|sed 's|"||g')
if [ "$license_server_location" == "null" ] || [ "$license_server_location" == "" ]; then
	#try global location
	set +e
	license_server_location=$(curl -s https://sentieon-bundle.s3.amazonaws.com/dnanexus/DNAnexus_app|jq '.license_server_location'|sed 's|"||g')
	set -e
	if [ "$license_server_location" == "null" ] || [ "$license_server_location" == "" ]; then
		license_server_location=master.sentieon.com:9010
	fi
fi
export SENTIEON_LICENSE=$license_server_location
export SENTIEON_AUTH_MECH=dnanexus_app
export SENTIEON_AUTH_DATA=~/Sentieon_License_encrypt
export SENTIEON_JOB_TAG=$job_tag

# Error checking: Check license validity before running anything
${SENTIEON_INSTALL_DIR}/bin/sentieon licclnt ping -s $SENTIEON_LICENSE 2> >(tee lic_errlog) || error_report "There is an issue with your license. Please contact Sentieon at support@sentieon.com and report your username, org, the billing org you are using to run Sentieon, and the following error message \"$(cat lic_errlog)\". The app will exit now."

#other settings
#export LD_PRELOAD=$SENTIEON_INSTALL_DIR/lib/libjemalloc.so.1
#export MALLOC_CONF=lg_dirty_mult:-1

```

</details>

The licence workflow proceeds as follows:

1. **`helper_funcs.sh:find_license_project`** — Identifies the correct DNAnexus project holding the Sentieon licence record. It searches for projects owned by `user-sentieon_license`, tagged `Sentieon_License`, first shared with the billing org, then with collaborator orgs, then with the launching user. Non-expired purchased licences are preferred over EVAL licences.

2. **`helper_funcs.sh:download_license_token`** — Downloads the licence record's JSON details via `dx describe`, base64-encodes them to `~/Sentieon_License_encrypt` (an atomic write via a temp file to avoid races).

3. **`license_auth.sh`** — Runs in the background as a refresh daemon (`while sleep 300; do download_license_token; done`), keeping the license token current throughout the 24-hour job window.

4. **`license_setup.sh`** — Resolves the licence server address (from the token, or from a Sentieon S3 endpoint, or falls back to `master.sentieon.com:9010`), then exports:
   - `SENTIEON_LICENSE` — host:port of the license server
   - `SENTIEON_AUTH_MECH=dnanexus_app` — authentication mechanism
   - `SENTIEON_AUTH_DATA` — path to the base64-encoded token
   - `SENTIEON_JOB_TAG` — job identifier for audit logging

5. **Licence ping** — `sentieon licclnt ping` is called before any compute starts; an informative error message is reported and the job fails fast if the licence is invalid.

## Bundled Software

All bioinformatics tools are shipped inside `resources/`, so no internet-based package installation is needed for the tools themselves.

<details>
<summary>Sentieon libexec binaries and shared libraries</summary>

```bash
ls resources/usr/local/sentieon-genomics-202503.01/libexec/
```

```output
bamslice
bwa
driver
fqidx
licclnt
licsrvr
LongReadUtil
minimap2
plot
pyexec
rcat
STAR
umi
util
```

```bash
ls resources/usr/local/sentieon-genomics-202503.01/lib/*.so* | xargs -I{} basename {}
```

```output
libblas.so.3
libblas.so.3.2
libblas.so.3.2.1
libbz2.so.1
libbz2.so.1.0.4
libgfortran.so.3
libgfortran.so.3.0.0
libgomp.so.1
libgomp.so.1.0.0
libhdf5.so.6
libhdf5.so.6.0.4
liblapack.so.3
liblapack.so.3.2
liblapack.so.3.2.1
liblzma.so.0
liblzma.so.0.0.0
libz.so.1
libz.so.1.2.8
```

</details>

| Component | Path | Role |
|---|---|---|
| `sentieon` (CLI wrapper) | `resources/usr/local/sentieon-genomics-202503.01/bin/sentieon` | Top-level entry point dispatching to libexec workers |
| `bwa` | `libexec/bwa` | Sentieon's optimised BWA-MEM implementation |
| `umi` | `libexec/umi` | UMI extraction worker |
| `driver` | `libexec/driver` | Sentieon's parallel driver (LocusCollector, Dedup, metrics algos) |
| `util` | `libexec/util` | BAM utilities including the sort used here |
| `licclnt` / `licsrvr` | `libexec/` | License client and server |
| `bamslice`, `fqidx`, `rcat`, `plot`, `pyexec` | `libexec/` | Supporting utilities (BAM slicing, FASTQ indexing, plotting, Python execution) |
| `minimap2`, `STAR`, `LongReadUtil` | `libexec/` | Long-read / RNA-seq aligners bundled with the suite but not used by this app |
| `samtools` | `resources/usr/bin/samtools` | Used exclusively to generate the `.fai` FASTA index |
| Shared libs | `lib/*.so` | BLAS, LAPACK, HDF5, zlib, bzip2, lzma, OpenMP, gfortran runtime (pinned versions to avoid system conflicts) |
| Python VCF lib | `lib/python/sentieon/vcflib/` | Python utilities for VCF manipulation (used by other Sentieon apps, present but inactive here) |

## Data Flow Summary

The complete end-to-end data flow for a paired-end, streaming run with two FASTQ pairs and metrics enabled:

```output
DNAnexus storage
  ├── R1_1.fq.gz, R2_1.fq.gz  (streamed via dx cat)
  ├── R1_2.fq.gz, R2_2.fq.gz  (streamed via dx cat)
  ├── genome.fa.gz             (streamed → genome/ref.fa + genome/ref.fa.fai)
  └── ref.bwa-index.tar.gz    (downloaded → genome/ref.fa.{amb,ann,bwt,pac,sa})

Worker local disk
  ├── [pair 0] sentieon umi extract | sentieon bwa mem | sentieon util sort → sorted_0.bam
  ├── [pair 1] sentieon umi extract | sentieon bwa mem | sentieon util sort → sorted_1.bam
  │
  ├── sentieon driver LocusCollector --consensus --umi_tag XR → score.txt.gz
  │     └── also: GCBias, MeanQualityByCycle, QualDistribution, InsertSize,
  │               AlignmentStat, CoverageMetrics  (if output_metrics=true)
  │
  ├── sentieon driver Dedup --score_info score.txt.gz → sorted.bam + sorted.bam.bai
  │     └── dedup_metrics.txt (always)
  │
  └── sentieon plot metrics → metrics.pdf  (if output_metrics=true)

DNAnexus storage (outputs)
  ├── <sample>_sorted.bam          (mappings_bam)
  ├── <sample>_sorted.bam.bai      (mappings_bam_bai)
  ├── <sample>_metrics/            (metrics — optional)
  │     ├── *.GCBias_metrics.txt
  │     ├── *.GCBiasSummary_metrics.txt
  │     ├── *.MeanQualityByCycle_metrics.txt
  │     ├── *.QualDistribution_metrics.txt
  │     ├── *.InsertSize_metrics.txt
  │     ├── *.AlignmentStat_metrics.txt
  │     ├── *.dedup_metrics.txt
  │     ├── *.coverage*
  │     └── *.metrics.pdf
  └── *.md5                        (md5sums — optional)
```

## Read Group Handling

Read group (`@RG`) tags are embedded in the output BAM at alignment time. If no `rg_info_csv` is provided, the app auto-generates sensible defaults:

| RG tag | Auto-generated value |
|--------|----------------------|
| `ID` | FASTQ filename prefix (stripped of `_R1`/`_1` suffixes) + lane number |
| `LB` | Sample name |
| `PU` | Same as `ID` |
| `PL` | Value of `read_group_platform` input (default: `ILLUMINA`) |
| `SM` | Value of `sample` input (or FASTQ-derived if not given) |

If `rg_info_csv` is supplied (a four-column CSV: filename, RG ID, RG LB, RG PU), the app validates that every R1 FASTQ filename has a corresponding row before proceeding, reporting a clear error if any are missing.

## Key Design Decisions

**Streaming I/O for speed**: By piping FASTQs and the genome FASTA rather than downloading them in full, the app overlaps network I/O with compute, substantially reducing wall-clock time on large samples.

**Single-pass metrics**: Collecting quality metrics in the same `LocusCollector` driver invocation as the consensus-score computation avoids an additional full scan of the BAMs.

**Configurable compression**: `bam_compression` (1/6/9) trades CPU time for intermediate storage. At the default of 1 (fastest), the 36-core instance prioritises throughput over disk efficiency.

**UMI tag propagation**: Using `sentieon bwa mem -C` ensures that custom SAM tags written by `umi extract` (including `XR` for the UMI sequence) flow through to the sorted BAM and are then read by `LocusCollector` via `--umi_tag XR`.

**Asynchronous upload with error resilience**: The upload machinery runs `dx upload` in the background while outputs are generated, and uses sentinel files rather than blocking to check completion, allowing multiple files to upload in parallel.

**Memory-bounded BWA**: The BWT in-memory size is capped at 400 GB (`bwt_max_mem`) with 6 GB headroom, preventing OOM on the high-memory 36-core instance while still loading the full human genome index into RAM.
