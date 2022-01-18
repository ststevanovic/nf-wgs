process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs", emit: fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

//
// Test with single-end data
//
// workflow test_fastqc_single_end {
//     input = [ [ id:'test', single_end:true ], // meta map
//               [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true) ]
//             ]

//     FASTQC ( input )
// }

// //
// // Test with paired-end data
// //
// workflow test_fastqc_paired_end {
//     input = [ [id: 'test', single_end: false], // meta map
//               [ file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
//                 file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true) ]
//             ]

//     FASTQC ( input )
// }