process HTML_COPYDIR {
   input:
    path sourceDir            // The input directory to be copied
    val  targetName           // The name of the output directory

    output:
    path "${targetName}"     // The result: a directory with a new name

    script:
    """
    mkdir -p ${targetName}
    cp -r ${sourceDir}/* ${targetName}/
    """
}
