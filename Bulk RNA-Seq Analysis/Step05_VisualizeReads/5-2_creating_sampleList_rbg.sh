#Go to directory with bigwig files
cd /tgc/TGCore_User_Data/WebData/cotney/hubs/gilmore/PWS_RNAseq
#Make a list of samples
ls > sample_list.txt
#Remove trailing part of sample names
sed -i 's/.str1.bw//g' sample_list.txt
sed -i 's/.str2.bw//g' sample_list.txt
#Remove duplicates
sort -u sample_list.txt > sample.txt
#Rename file
mv sample.txt sample_list.txt
