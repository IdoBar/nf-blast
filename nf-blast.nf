#!/usr/bin/env nextflow

nextflow.enable.dsl = 2 
/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 */
 

def helpMessage() {

     log.info """
      Usage:
      The typical command for running the pipeline is as follows:
      nextflow run blast-nf.nf --app blastn --query QUERY.fasta --db "blastDatabaseDirectory/blastPrefixName" -profile local
      nextflow run blast-nf.nf --app "diamond blastp" --query QUERY.faa --db "diamondDirectory/diamondDB.dmnd" -profile local

      Mandatory arguments:
       --app                          BLAST/DIAMOND program to use (diamond blastp/x must be quoted!)
                                      Valid options: [blastn, blastp, tblastn, blastx, 'diamond blastp', 'diamond blastx'] 
       --query                        Query fasta file of sequences you wish to BLAST
       --db                           Prefix name of the BLAST or DIAMOND database (full path required). 
                                      If blast database is provided for diamond and taxonomy information is requested
                                      then a diamond database will be created (see Taxonomy options below). 
                                      Default: [$BLASTDB env variable] 
       -profile                       Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, apptainer, singularity, local, blastn_tax, diamond_tax, test, test_tax,
                                      test_p, test_p_tax, test_d, test_d_tax

       Optional arguments:
       --out                          Output filename of final BLAST output. Default: [QUERY.outfmt6]
       --outDir                       Output folder for the results. Default: [blast-nf-results]
       --outfmt                       Output format (must be quoted!). Default: ['6 std']
       --blastOpts                    Additional options for BLAST command (must be quoted!). Default: ['-evalue 1e-10 -max_target_seqs 20']
       --dmndOpts                     Additional options for BLAST command (must be quoted!). Default: ['-e 1e-10 -k 20'] 
       --chunkSize                    Number of fasta records to use in each job when splitting the query fasta file. Default: [500]
                                      This option can also take the size of each subquery (like 200.KB, 5.KB, etc.) 
       --queueSize                    Maximum number of jobs to be queued [20]

       Taxonomy options:
       --taxDbDir                     Location of taxonomy db files (prot.accession2taxid.FULL.gz, nodes.dmp and names.dmp) to allow Diamond 
                                      to return taxonomic information columns. If the required files cannot be found in the path 
                                      they will be automatically downloaded from the web.
                                      Information about the required files and where to download them can be found at 
                                      https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#makedb-options
                                      Default: [same path as the database]
        --taxListFile                 A file with list of taxonomy IDs to limit the search space.
       --help                         This usage statement.
     """
 }



 // Show help message
 if (params.help) {
     helpMessage()
     exit 0
 }

def db_name = file(params.db).name
// diamond_db = file(params.db).name
// db_path = (params.app =~ /D$/ && params.outfmtString =~ /staxids|sscinames|scomnames|sskingdoms/ && file("${params.db}.dmnd").exists()) ? "${params.db}.dmnd" : file(params.db).name
def db_basename = file(params.db).baseName
def db_dir = file(params.db).parent
def db_prefix = db_dir.resolve("${db_basename}")
def tax_db_dir = params.taxDbDir ?: db_dir
// out_dir = "${params.outDir}/${params.app}"
// Check if the chunks are provided as Memory units and if not assume KB
def chunk_size = params.chunkSize // =~ /\d+\.*\w[bB]$/ ? MemoryUnit.of( "${params.chunkSize}.KB" ) : MemoryUnit.of( params.chunkSize )
if (params.chunkSize ==~ /\d+\.*\w[bB]$/) {
    (mem_value, mem_suffix) = (chunk_size =~ /(\d+)\.*(\w[bB])$/)[0]
    chunk_size = mem_value + "." + mem_suffix.toUpperCase()
}

// Format output filename with as "query"
def output_fmt = params.outfmtString =~ /(\d+) (.+)/ 
def tax_filt = file(params.taxListFile)
def app_str = params.app.replaceAll(/\s/, "-")
def res_name = params.out ?: "${file(params.query).baseName}.${app_str}.${file(params.db).simpleName}.outfmt${output_fmt[0][1]}"
def db_ready = false
// dag_file = (nextflow.version.matches('22.04+')) ? "${params.outDir}/${params.tracedir}/blast-nf_dag.mmd" : "${params.outDir}/${params.tracedir}/blast-nf_dag.html"
/*
PROCESSES
*/
process blast {
    label 'blast'

    publishDir params.outDir, mode: 'copy', pattern: 'software_versions.txt'

    input:
        path 'query.fa'
        val dbPrefix 
       // val dbName
        path taxlist
 
    output:
        path 'blast_results'
       // path 'software_versions.txt'

    script:
    def taxid_filt = taxlist.name != 'NO_FILE' ? "-taxidlist ${taxlist}" : ''
    """
    ${params.app} -db ${dbPrefix} -query 'query.fa' -outfmt '${params.outfmtString}' ${taxid_filt} -num_threads ${task.cpus} ${params.blastOpts} > blast_results
    blastn -version | head -n1 > software_versions.txt
    """
}

process taxids_file2list {
    executor 'local'
    input:
        path taxid_file

    output:
        stdout

    shell:
    '''
    tax_str=!( gawk -v ORS="," '1' !{taxid_file} | sed 's/,$//' )
    echo "--taxonlist ${tax_str}"
    '''
}


process extract_db_fasta {
    label 'blast'
    label 'sc_medium'
    input:
        val dbPrefix 

    output:
        path "db.fasta" 

    """
    blastdbcmd -entry all -db ${dbPrefix} -out db.fasta
    """

}

process download_tax_db {
    label 'download'

    publishDir tax_db_dir, mode = 'copy' 

    output:
        path "*.dmp"
        path "prot.accession2taxid.FULL.gz"

    script:
    """
    aria2c -x5 -c --auto-file-renaming=false https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
    aria2c -x5 -c --auto-file-renaming=false https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    unzip taxdmp.zip
    """
}

process diamond_prep_tax_db {
    label 'diamond'

    // publishDir "${tax_db_dir}", mode = 'copy' 

    input:
        val dbPrefix
        path fasta
        path acc2tax
        path names
        path nodes
    
    output:
        val true  
    
    script:
    """
    diamond makedb -p ${task.cpus} --in ${fasta} -d ${dbPrefix} --taxonmap ${acc2tax} --taxonnodes ${nodes} --taxonnames ${names}
    """

}

process diamond_prep_blast_db {
    label 'diamond'
    // publishDir "${dbDir}", mode = 'copy' 

    input:
        // val dbDir
        // val dbName
        val dbPrefix // from "${db_name}"
    
    output:
        val true  
    
    script:
    """
    diamond prepdb -p ${task.cpus} -d ${dbPrefix} 
    """

}


process diamond {
    label 'diamond'

    publishDir params.outDir, mode: 'copy', pattern: 'software_versions.txt'

    input:
        path 'query.fa' // from ch_fasta
        val dbPrefix //from diamond_db_ch
        // val dbName
        val tax_filt_string
        val ready
    output:
        path 'diamond_results'
       // path 'software_versions.txt'

    script:
    """
    ${params.app} -d ${dbPrefix} -q 'query.fa' ${tax_filt_string} -f ${params.outfmtString} -p ${task.cpus} ${params.dmndOpts} -o diamond_results
    diamond version > software_versions.txt
    """
} 

/*

workflow diamond_prep_flow { 
    take: 
        fmt_string
        db_prefix
        tax_db_dir
        tax_filt
        // db_ready
    main:
        if (fmt_string =~ /staxids|sscinames|sphylums|skingdoms/) {
            // def diamond_db = db_dir.resolve("${db_basename}.dmnd")
            if (file("${db_prefix}.dmnd").isEmpty()) {
                def acc2taxid = file("${tax_db_dir}/prot.accession2taxid.FULL.gz")
                def names_dmp = file("${tax_db_dir}/names.dmp")
                def nodes_dmp = file("${tax_db_dir}/nodes.dmp")
                db_fasta_ch = file("${db_prefix}.fasta").isEmpty() ? extract_db_fasta(db_prefix) : Channel.fromPath("${db_prefix}.fasta")
                if (acc2taxid.isEmpty() || names_dmp.isEmpty() || nodes_dmp.isEmpty()) 
                    download_tax_db()
                
                diamond_prep_tax_db(db_prefix, db_fasta_ch, Channel.fromPath(acc2taxid), Channel.fromPath(names_dmp), Channel.fromPath(nodes_dmp))
            } 
        } 
        else if (file("${db_prefix}.acc").isEmpty()) 
                diamond_prep_blast_db(db_prefix)
        // process taxonomy filtering string
        def tax_filt_d (tax_filt_file) {
            if (tax_filt_file.name == 'NO_FILE')
                return ''
            tax_file_str = tax_filt_file.text
            tax_str = tax_file_str.replaceAll(/\n/, ",")
            tax_str = "--taxonlist " + tax_str.replaceFirst(/,$/, "")
            return tax_str
            //tax_filt_file.text.replaceAll(/\n/, ",").replaceFirst(/,$/, "")        
        }
        tax_filt_string = tax_filt_d(tax_filt)
        // tax_filt.name != 'NO_FILE' ? "--taxonlist ${taxids_file2list(tax_filt)}" : ''
            
        db_ready = true
        
    emit:
        tax_filt_string
        db_ready

}
*/
workflow {
    
    // Define your channels
    if (chunk_size ==~ /\d+\.\w+[bB]/) {
        Channel.fromPath(params.query)
                .splitFasta(size: MemoryUnit.of( chunk_size ), file: true)
                .set { ch_fasta }
    }
    else {
        Channel.fromPath(params.query)
                .splitFasta(by: chunk_size, file: true)
                .set { ch_fasta }
    }
    
    
    
    // db_ch = Channel.fromPath(params.db)
    // Run the homology search process 
    if (params.app =~ /diamond/) {
        dmnd_db = "${db_prefix}"
        tax_str_ch = tax_filt.name == 'NO_FILE' ? '' : taxids_file2list(tax_filt) 
        if (params.outfmtString =~ /staxids|sscinames|sphylums|skingdoms/) {
            // def diamond_db = db_dir.resolve("${db_basename}.dmnd")
            if (file("${db_prefix}.dmnd").exists())
                db_ready_ch = channel.value(true)
            else {
                acc2taxid = file("${tax_db_dir}/prot.accession2taxid.FULL.gz")
                names_dmp = file("${tax_db_dir}/names.dmp")
                nodes_dmp = file("${tax_db_dir}/nodes.dmp")
                // download acc2tax if missing
                if (acc2taxid.isEmpty() || names_dmp.isEmpty() || nodes_dmp.isEmpty()) 
                    download_tax_db()
                db_fasta_ch = file("${db_prefix}.fasta").isEmpty() ? extract_db_fasta(db_prefix) : Channel.fromPath("${db_prefix}.fasta")
                diamond_prep_tax_db(db_prefix, db_fasta_ch, acc2taxid, names_dmp, nodes_dmp)
                db_ready_ch = diamond_prep_tax_db.out
            }
            dmnd_db = "${db_prefix}.dmnd"
        } 
        else if (file("${db_prefix}.acc").exists()) 
            db_ready_ch = channel.value(true)
        else {
            diamond_prep_blast_db(db_prefix)
            db_ready_ch = diamond_prep_blast_db.out
        } 
     
        ch_hits = diamond(ch_fasta, dmnd_db, tax_str_ch , db_ready_ch)      
        // tax_filt_ch = Channel.fromPath(params.taxListFile)        
        // run diamond
        // (tax_filt_ch, db_ready_ch) = diamond_prep_flow(Channel.value(db_prefix), Channel.value(tax_db_dir), Channel.value(tax_filt), Channel.value(db_ready))
        // if (diamond_prep_flow.out[1]) {
        
        
        //}
            
        /*
        if (diamond_prep_flow.out[1])
            ch_hits = diamond(ch_fasta, db_prefix, diamond_prep_flow.out[0], diamond_prep_flow.out[1])
            */
    }
    else {
        ch_hits = blast(ch_fasta, db_prefix, tax_filt)
    }
         
          
    // Collect results
    ch_hits
        .collectFile(name: res_name, storeDir: params.outDir)
        .subscribe { println "Entries are saved to file: $it" }
}

