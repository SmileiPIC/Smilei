#!/bin/sh

outdir="."
command="smilei"

SMILEI_DIR=$( cd "$( dirname "$0" )" && pwd )

usage() { echo "Usage: $0 [-d output directory] command namelist" 1>&2; exit 1; }

while getopts "hd:" o; do
    case "${o}" in
        d)
            outdir=${OPTARG}
            ;;
        c)
            command=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ $# -eq 1 ]
  then
    mkdir -p "${outdir}"
    git --git-dir=${SMILEI_DIR}/.git archive --format=tar.gz HEAD > "${outdir}/code.tgz"
    cd $outdir
    $command $1
    echo "to be continued"
  else
    usage
fi



