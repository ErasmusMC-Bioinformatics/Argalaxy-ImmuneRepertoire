set -e

dir="$(cd "$(dirname "$0")" && pwd)"

input=$1
species=$2
locus=$3
output=$4

declare -A speciesdict

speciesdict=(["Rattus norvegicus functional"]="rat" ["Rattus norvegicus non-functional"]="rat" ["Oryctolagus cuniculus functional"]="rabbit" ["Oryctolagus cuniculus non-functional"]="rabbit" ["Mus musculus functional"]="mouse" ["Mus musculus non-functional"]="mouse" ["Homo sapiens functional"]="human" ["Homo sapiens non-functional"]="human" ["Macaca mulatta non-functional"]="rhesus_monkey" ["Macaca mulatta functional"]="rhesus_monkey")

echo "Species: $species ${speciesdict[$species]}"

species="${speciesdict[$species]}"

if [ "$species" == "" ]
then
	>&2 echo "Species not possible with igBLASTn, use IMGT"
	exit 1
fi

echo "$input $species $locus $output"

java -Xmx16G -jar $IGBLASTWRP/igblastwrp.jar -p 4 -S $species -R $locus ${input} $PWD/blasted_output 2>&1

Rscript --verbose $dir/igblast.r "$PWD/blasted_output.L2.txt" "$output" 2>&1
