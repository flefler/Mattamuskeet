# Lets make some quick trees for our taxa

## Packages used 
[IQ-TREE ](https://github.com/iqtree/iqtree2)

[GToTree](https://github.com/AstrobioMike/GToTree)

# Gather the genomes from our internal database
```
conda activate gtotree
cd /blue/hlaughinghouse/flefler/MATAMUSK

mkdir GENOMES
mkdir LOGS

MAGs=`grep 'c__Cyanobacteriia' /orange/hlaughinghouse/flefler/BLCC_Genomes/08_gtdbtk/gtdbtk.bac120.summary.tsv | grep 'F208|F209|F210|F211|F215|F218' | cut -f 1`
for MAG in $MAGs; do
  cp /orange/hlaughinghouse/flefler/BLCC_Genomes/05_GENOMES/${MAG}.fa.gz GENOMES/${MAG}.fa.gz
done
```
## List genomes in file
```
mkdir MICRO 
mkdir UMA
mkdir PELTO
mkdir RAPHIDO
echo "F209_bin4.1_re.2_re.fa.gz F210_bin.1.polca.fa.gz" > MICRO/micro_list.txt
echo "F208_bin6.fa.gz F215_bin5.fa.gz" > UMA/uma_list.txt
echo "F211_bin3.1_re.fa.gz" > PELTO/Pleto_list.txt
echo "F218_SemiBin_2.fa.gz" > RAPHIDO/raphidio_list.txt
```

# Get genomes from GTDB/NCBI
probably could have done a loopy-de-loop here, but alas... no
```
cd MICRO && gtt-get-accessions-from-GTDB -t Microcystis --GTDB-representatives-only #There is a lot so lets refine this
cd UMA && gtt-get-accessions-from-GTDB -t Umezakia 
cd PELTO && gtt-get-accessions-from-GTDB -t Pelatocladus &&\
            gtt-get-accessions-from-GTDB -t Fischerella --GTDB-representatives-only &&\
            cat *.tsv > combined.tsv && cat *.txt > combined.txt
cd RAPHIDO && gtt-get-accessions-from-GTDB -t Raphidiopsis
```
# TREE TIME, almost

## Identify single copy genes
```
N=GTOTREE_MICRO
CMD="cd MICRO && GToTree -f micro_list.txt -a GTDB-Microcystis-genus-GTDB-rep-accs.txt -H Cyanobacteria -B -D -n 2 -j 8 -N"
sbatch -A hlaughinghouse -J ${N} -c 8 --mem=10G -o LOGS/${N}_o.txt -e LOGS/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"

N=GTOTREE_UMA
CMD="cd UMA && GToTree -f uma_list.txt -a GTDB-Umezakia-genus-accs.txt -H Cyanobacteria -B -D -n 2 -j 8 -N"
sbatch -A hlaughinghouse -J ${N} -c 8 --mem=10G -o LOGS/${N}_o.txt -e LOGS/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"

N=GTOTREE_PELTO
CMD="cd PELTO && GToTree -f Pleto_list.txt -a combined.txt -H Cyanobacteria -B -D -n 2 -j 8 -N"
sbatch -A hlaughinghouse -J ${N} -c 8 --mem=10G -o LOGS/${N}_o.txt -e LOGS/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"

N=GTOTREE_RAPHIDO
CMD="cd RAPHIDO && GToTree -f raphidio_list.txt -a GTDB-Raphidiopsis-genus-accs.txt -H Cyanobacteria -B -D -n 2 -j 8 -N"
sbatch -A hlaughinghouse -J ${N} -c 8 --mem=10G -o LOGS/${N}_o.txt -e LOGS/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"
```
## Tree, for real this time
```
N=GTOTREE_MICRO
CMD="cd MICRO && iqtree2 -s GToTree_output/Aligned_SCGs_mod_names.faa --seed 42069 --mem 100G -T 16 -m MFP --merit AICc -B 1000 -pre MICRO"
sbatch -A hlaughinghouse -J ${N} -c 16 --mem=100G -o LOGS/${N}_o.txt -e LOGS/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"

N=GTOTREE_UMA
CMD="cd UMA && iqtree2 -s GToTree_output/Aligned_SCGs_mod_names.faa --seed 42069 --mem 100G -T 16 -m MFP --merit AICc -B 1000 -pre UMA"
sbatch -A hlaughinghouse -J ${N} -c 16 --mem=100G -o LOGS/${N}_o.txt -e LOGS/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"

N=GTOTREE_PELTO
CMD="cd PELTO && iqtree2 -s GToTree_output/Aligned_SCGs_mod_names.faa --seed 42069 --mem 100G -T 16 -m MFP --merit AICc -B 1000 -pre PELTO"
sbatch -A hlaughinghouse -J ${N} -c 16 --mem=100G -o LOGS/${N}_o.txt -e LOGS/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 3:00:00 --wrap="${CMD}"

N=GTOTREE_RAPHIDO
CMD="cd RAPHIDO && iqtree2 -s GToTree_output/Aligned_SCGs_mod_names.faa --seed 42069 --mem 100G -T 16 -m MFP --merit AICc -B 1000 -pre RAPHIDO"
sbatch -A hlaughinghouse -J ${N} -c 16 --mem=100G -o LOGS/${N}_o.txt -e LOGS/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=flefler@ufl.edu -t 1:00:00 --wrap="${CMD}"
```
