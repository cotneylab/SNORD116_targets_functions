cat sample_list.txt | awk '{ \
print "\ntrack type=bigWig name=%"$1"_R% description=%PWS_RNAseq_R% bigDataUrl=http://URL/"$1".str1.bw color=255,0,0 visibility=full yLineOnOff=on autoScale=on yLineMark=0 alwaysZero=off viewLimits=2:15 graph Type=bar maxHeightPixels=64:32:16 windowingFunction=mean+whiskers smoothingWindow=off negateValues=on" \
}' > UCSC_upload.txt
cat sample_list.txt | awk '{ \
print "\ntrack type=bigWig name=%"$1"_F% description=%PWS_RNAseq_F% bigDataUrl=http://URL/"$1".str2.bw color=0,0,255 visibility=full yLineOnOff=on autoScale=on yLineMark=0 alwaysZero=off viewLimits=2:15 graph Type=bar maxHeightPixels=64:32:16 windowingFunction=mean+whiskers smoothingWindow=off negateValues=off" \
}' >> UCSC_upload.txt
sed -i 's/%/"/g' UCSC_upload.txt
