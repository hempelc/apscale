import subprocess, datetime, gzip, os, pickle, glob, openpyxl, shutil, psutil, re, sys
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from io import StringIO
from tqdm import tqdm
from openpyxl.utils.dataframe import dataframe_to_rows


## denoising function to denoise all sequences the input fasta with a given alpha and minsize
def denoise(
    project=None, comp_lvl=None, cores=None, alpha=None, minsize=None, coi=None
):
    """Function to apply denoisind to a given gzipped file. Outputs a fasta file with all
    centroid sequences."""

    ## define the name for the output fasta
    ## create an output path to write to
    sample_name_out_1 = "ESVs_with_chimeras.fasta.gz"
    gz_fasta = Path(project).joinpath(
        "6_dereplication_pooling",
        "data",
        "pooling",
        "pooled_sequences_dereplicated.fasta.gz",
    )
    fasta = os.path.splitext(gz_fasta)[0]
    fasta_edited = f"{os.path.splitext(fasta)[0]}_seqidedit.fasta"
    output_path = Path(project).joinpath("8_denoising", "data", sample_name_out_1)
    temp_path = Path(project).joinpath("8_denoising", "temp")

    # Open the gzip-compressed file and the output file
    with gzip.open(gz_fasta, "rb") as gz_file, open(fasta, "wb") as output:
        # Read the compressed data and write it to the output file
        output.write(gz_file.read())

    # collect number of processed reads
    seqs = len([1 for line in open(fasta) if line.startswith(">")])

    if coi == True:
        ## give user output
        print(
            f'{datetime.datetime.now().strftime("%H:%M:%S")}: Starting denoising with DnoisE. This may take a while.'
        )

        # Replace sequence IDs with ascending numbers to remove duplicate IDs (otherwise DnoisE fails)
        sequence_number = 1
        with open(fasta, "r") as infile, open(fasta_edited, "w") as outfile:
            for line in infile:
                if line.startswith(">seq:"):
                    # Extract the size information
                    size_info = line.split(";")[1]
                    # Replace the sequence number in the ID and keep the size information
                    line = f">seq:{sequence_number};{size_info}"
                    # Increment the sequence number
                    sequence_number += 1
                outfile.write(line)

        # run DnoisE to denoise reads, remove unused denoising info, and rename DnoisE fasta
        # min_abun is set to 8 (minimum read abundance) to match the default unoise setting
        with open(temp_path.joinpath("dnoise_log.txt"), "w") as output:
            f = subprocess.run(
                [
                    "dnoise",
                    "--fasta_input",
                    fasta_edited,
                    "--fasta_output",
                    output_path.with_suffix(""),
                    "--min_abun",
                    str(minsize),
                    "-y",
                    "--cores",
                    str(cores),
                ],
                stdout=output,
            )
        os.remove(f'{output_path.with_suffix("")}_Adcorr_denoising_info.csv')
        os.rename(
            f'{output_path.with_suffix("")}_Adcorr_denoised_ratio_d.fasta',
            output_path.with_suffix(""),
        )

        # Remove uncompressed fasta file
        os.remove(fasta_edited)

    elif coi == False:
        ## give user output
        print(
            f'{datetime.datetime.now().strftime("%H:%M:%S")}: Starting denoising with vsearch unoise. This may take a while.'
        )

        with open(output_path.with_suffix(""), "w") as output:
            f = subprocess.run(
                [
                    "vsearch",
                    "--cluster_unoise",
                    fasta,
                    "--unoise_alpha",
                    str(alpha),
                    "--minsize",
                    str(minsize),
                    "--sizein",
                    "--sizeout",
                    "--centroids",
                    "-",
                    "--fasta_width",
                    str(0),
                    "--quiet",
                    "--log",
                    temp_path.joinpath("denoising_log.txt"),
                    "--threads",
                    str(cores),
                ],
                stdout=output,
                stderr=subprocess.DEVNULL,
            )

    else:
        print(
            f'"coi" must be set to either True or False, current setting: {coi}',
            file=sys.stderr,
        )
        sys.exit()

    # Collect number of ESVs
    esvs = len(
        [1 for line in open(output_path.with_suffix("")) if line.startswith(">")]
    )

    ## compress the output, remove uncompressed output
    with open(output_path.with_suffix(""), "rb") as in_stream, gzip.open(
        output_path, "wb", comp_lvl
    ) as out_stream:
        shutil.copyfileobj(in_stream, out_stream)
    os.remove(output_path.with_suffix(""))

    print(
        "{}: Denoised unique {} sequences into {} ESVs.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), seqs, esvs
        )
    )
    print(
        "{}: Starting chimera removal from the ESVs. This may take a while.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )

    ## run vsearch --uchime_denovo to remove chimeric sequences from the ESVs
    f = subprocess.run(
        [
            "vsearch",
            "--uchime_denovo",
            Path(project).joinpath("8_denoising", "data", sample_name_out_1),
            "--relabel",
            "ESV_",
            "--nonchimeras",
            Path(project).joinpath(
                "8_denoising", "{}_ESVs.fasta".format(Path(project).stem)
            ),
            "-fasta_width",
            str(0),
            "--quiet",
        ]
    )

    ## collect processed and passed reads from the output fasta, since it is not reported in the log
    f = list(
        SimpleFastaParser(
            open(
                Path(project).joinpath(
                    "8_denoising", "{}_ESVs.fasta".format(Path(project).stem)
                )
            )
        )
    )
    print(
        "{}: {} chimeras removed from {} ESV sequences.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), int(esvs) - len(f), esvs
        )
    )
    print(
        "{}: ESVs saved to {}.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            Path(project).joinpath(
                "7_otu_clustering", "{}_OTUs.fasta".format(Path(project).stem)
            ),
        )
    )


## remapping function to remap the individual reads to the ESVs via vsearch
def remapping_esv(file, project=None):
    """Function to remap the sequences of a dereplicated file against the ESV list
    as database."""

    ## extract the sample name from the file name for the otu table
    sample_name_out = "{}".format(
        Path(file).with_suffix("").with_suffix("").name
    ).replace("_PE_trimmed_filtered_dereplicated", "")

    ## run vsearch --search_exact to remap the individual files vs the generated
    ## ESV fasta, capture log and directly pickle the output as dataframe for read table generation
    f = subprocess.run(
        [
            "vsearch",
            "--search_exact",
            Path(file),
            "--db",
            Path(project).joinpath(
                "8_denoising", "{}_ESVs.fasta".format(Path(project).stem)
            ),
            "--output_no_hits",
            "--maxhits",
            "1",
            "--otutabout",
            "-",
            "--quiet",
            "--threads",
            str(1),
            "--log",
            Path(project).joinpath(
                "8_denoising", "temp", "{}_mapping_log.txt".format(sample_name_out)
            ),
        ],
        capture_output=True,
    )

    ## directly parse the output to a pandas dataframe
    esv_tab = pd.read_csv(StringIO(f.stdout.decode("ascii", errors="ignore")), sep="\t")

    ## handle empty outputs correctly
    if not esv_tab.empty:
        esv_tab = esv_tab.set_axis(["ID", sample_name_out], axis=1, copy=False)
    else:
        esv_tab[sample_name_out] = ""

    ## collect number of esvs from the output, pickle to logs
    with open(
        Path(project).joinpath(
            "8_denoising", "temp", "{}_mapping_log.txt".format(sample_name_out)
        )
    ) as log_file:
        content = log_file.read()
        exact_matches = re.findall("Matching query sequences: (\d+)", content)[0]
        version = re.findall("vsearch ([\w\.]*)", content)[0]
        finished = "{}".format(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))

    ## give user output
    print(
        "{}: {}: {} exact matches found ({} reads).".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            sample_name_out,
            exact_matches,
            esv_tab[sample_name_out].sum(),
        )
    )

    ## pickle log data first for log generation
    with open(
        Path(project).joinpath(
            "8_denoising", "temp", "{}_log.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump(
            [
                sample_name_out,
                finished,
                version,
                exact_matches,
                esv_tab[sample_name_out].sum(),
            ],
            log,
        )

    ## pickle otu tab dataframes for otu table generation
    with open(
        Path(project).joinpath(
            "8_denoising", "temp", "{}_esv_tab.pkl".format(sample_name_out)
        ),
        "wb",
    ) as log:
        pickle.dump(esv_tab, log)


## main function for the denoising script
def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Will denoise the dataset, perform chimera removal, remap the individual files and
    generate an ESV table."""

    ## create temporal output folder
    try:
        os.mkdir(Path(project).joinpath("8_denoising", "temp"))
    except FileExistsError:
        pass

    ## collect variables from the settings file
    gen_settings = pd.read_excel(
        Path(project).joinpath("Settings.xlsx"), sheet_name="0_general_settings"
    )
    cores, comp_lvl = (
        gen_settings["cores to use"].item(),
        gen_settings["compression level"].item(),
    )

    settings = pd.read_excel(
        Path(project).joinpath("Settings.xlsx"), sheet_name="8_denoising"
    )
    alpha, minsize, coi, to_excel = (
        settings["alpha"].item(),
        settings["minsize_denoising"].item(),
        settings["coi"].item(),
        settings["to excel"].item(),
    )

    ## denoise the dataset
    denoise(
        project=project,
        comp_lvl=comp_lvl,
        cores=cores,
        alpha=alpha,
        minsize=minsize,
        coi=coi,
    )

    ## gather files for remapping of ESVs
    input = glob.glob(
        str(
            Path(project).joinpath(
                "6_dereplication_pooling", "data", "dereplication", "*.fasta.gz"
            )
        )
    )

    print(
        "{}: Starting to remap {} input files.".format(
            datetime.datetime.now().strftime("%H:%M:%S"), len(input)
        )
    )

    ## run remapping parallelized to speed up the process
    Parallel(n_jobs=cores)(
        delayed(remapping_esv)(file, project=project) for file in input
    )

    ## write log for the denoising from pkl logs
    summary_logs = glob.glob(
        str(Path(project).joinpath("8_denoising", "temp", "*_log.pkl"))
    )
    summary = [pickle.load(open(line, "rb")) for line in summary_logs]

    log_df = pd.DataFrame(
        summary,
        columns=[
            "File",
            "finished at",
            "program version",
            "exact matches",
            "reads matched",
        ],
    )
    log_df = log_df.sort_values(by="File")
    log_df.to_excel(
        Path(project).joinpath("8_denoising", "Logfile_8_denoising.xlsx"),
        index=False,
        sheet_name="8_denoising",
    )

    ## add log to the project report
    with pd.ExcelWriter(
        Path(project).joinpath("Project_report.xlsx"),
        mode="a",
        if_sheet_exists="replace",
        engine="openpyxl",
    ) as writer:
        log_df.to_excel(writer, sheet_name="8_denoising", index=False)

    ## generate OTU table, first extract all OTUs and sequences from fasta file
    esv_list = list(
        SimpleFastaParser(
            open(
                Path(project).joinpath(
                    "8_denoising", "{}_ESVs.fasta".format(Path(project).stem)
                )
            )
        )
    )
    esv_table = pd.DataFrame(esv_list, columns=["ID", "Seq"])
    seq_col = esv_table.pop("Seq")

    ## extract individual ESV tabs from the clustering output, rename columns correctly, merge individual tabs
    esv_tabs = glob.glob(
        str(Path(project).joinpath("8_denoising", "temp", "*_esv_tab.pkl"))
    )
    esv_tabs = [pickle.load(open(tab_file, "rb")) for tab_file in esv_tabs]
    esv_tabs = [tab.rename(columns={tab.columns[0]: "ID"}) for tab in esv_tabs]
    esv_tabs = [
        pd.merge(esv_table, tab, on="ID", how="outer").set_index("ID")
        for tab in tqdm(esv_tabs, desc="Generating ESV table")
    ]

    ## collapse all individual dataframes into the ESV table, replace nan values with 0
    esv_table = pd.concat(esv_tabs, axis=1)
    esv_table = esv_table.reset_index(level=0).fillna(0)
    esv_table = pd.concat(
        [
            esv_table[["ID"]],
            esv_table[esv_table.columns.difference(["ID"])].sort_index(axis=1),
        ],
        ignore_index=False,
        axis=1,
    )

    ## move sequences to the end of the dataframe
    esv_table.insert(len(esv_table.columns), "Seq", seq_col)

    ## save the final OTU table if selected
    if to_excel:
        wb = openpyxl.Workbook(write_only=True)
        ws = wb.create_sheet("ESV table")

        ## save the output line by line for optimized memory usage
        for row in tqdm(
            dataframe_to_rows(esv_table, index=False, header=True),
            total=len(esv_table.index),
            desc="{}: Lines written to ESV table".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            ),
            unit=" lines",
        ):
            ws.append(row)

        ## save the output (otu table)
        print(
            "{}: Saving the ESV table to excel. This may take a while.".format(
                datetime.datetime.now().strftime("%H:%M:%S")
            )
        )
        wb.save(
            Path(project).joinpath(
                "8_denoising", "{}_ESV_table.xlsx".format(Path(project).stem)
            )
        )
        wb.close()
        print(
            "{}: ESV table saved to {}.".format(
                datetime.datetime.now().strftime("%H:%M:%S"),
                Path(project).joinpath(
                    "8_denoising", "{}_ESV_table.xlsx".format(Path(project).stem)
                ),
            )
        )

    print(
        "{}: Saving the ESV table to parquet. This may take a while.".format(
            datetime.datetime.now().strftime("%H:%M:%S")
        )
    )
    esv_table.to_parquet(
        Path(project).joinpath(
            "8_denoising", "{}_ESV_table.parquet.snappy".format(Path(project).stem)
        ),
        index=False,
    )
    print(
        "{}: ESV table saved to {}.".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            Path(project).joinpath(
                "8_denoising", "{}_ESV_table.parquet.snappy".format(Path(project).stem)
            ),
        )
    )

    ## remove temporary files
    shutil.rmtree(Path(project).joinpath("8_denoising", "temp"))


if __name__ == "__main__":
    main()
