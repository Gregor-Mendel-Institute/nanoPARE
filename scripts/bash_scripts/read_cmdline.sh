#!/bin/bash
############################
# READING THE COMMAND LINE #
############################

if [ $# -eq 0 ]; then
    echo "No commandline arguments parsed."
else
while [ "$1" != "" ]; do
    case $1 in
    -h | --help )               usage; exit 1
                                ;;
    --lmod )                    LMOD=1
                                ;;
    --ram )                     shift; RAM=$1
                                ;;
    --cpus )                    shift; CPUS=$1
                                ;;
    -G | --genome )             shift; GENOME_FASTA=$1
                                ;;
    -A | --annotation )         shift; ANNOTATION_GFF=$1
                                ;;
    -R | --reference )          shift; REFERENCE_TABLE=$1
                                ;;
    -L | --line )               shift; JOB_NUMBER=$1
                                ;;
    -N | --name )               shift; SAMPLE_NAME=$1;
                                ;;
    -T | --type )               shift; SAMPLE_TYPE=$1;
                                ;;
    --icomp )                   shift; ICOMP=$1
                                ;;
    --uug )                     shift; UUG=$1
                                ;;
    --upstream )                shift; UPSTREAM=$1
                                ;;
    --rpm )                     shift; RPM=$1
                                ;;
    --fraglen )                 shift; FRAGLEN=$1
                                ;;
    --bandwidth )               shift; BANDWIDTH=$1
                                ;;
    --kernel )                  shift; KERNEL=$1
                                ;;
    --mask )                    shift; MASK_NAME=$1
                                ;;
    * )                         echo "Argument not recognized."; usage; exit 1
                                ;;
    esac
    shift
done
fi

