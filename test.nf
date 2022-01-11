Channel
    .from("${params.outdir}/BWA/bwa_index/")
    .collectFile(name: 'test.txt', newLine: true)
    .subscribe {
        println "Entries are saved to file: $it"
        println "File content is: ${it.text}"
    }