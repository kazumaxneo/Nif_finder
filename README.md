# Nif_finder

#example code  
ls *faa > list  
cat list |while read line  
do  
python Nif_finderv0.13.py --query_file $line --output_prefix ${line%.faa}_result --cpu 8  \
--reference_files nifH/nifHclassification nifD/nifDclassification nifK/nifKclassification nifE/nifEclassification nifN/nifNclassification nifB/nifBclassification  \
 --target_files nifH/proteins_hmm nifD/proteins_hmm nifK/proteins_hmm nifE/proteins_hmm nifN/proteins_hmm nifB/proteins_hmm  
done
