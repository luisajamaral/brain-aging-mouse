#conda activate homer 
#for i in * ; do findMotifsGenome.pl $i mm10 motifs/${i%.bed} -nomotif -bg ../../../../../data/snATAC/peaks/9mo+18mo_summits_ext1k.bed ; done

#for i in * ; do findMotifsGenome.pl $i mm10 motifs/${i%.bed} -nomotif -bg ../../../../../data/snATAC/peaks/9mo+18mo_summits_ext1k.bed ; done

for i in grouped_peaks/*Blue* ; do findMotifsGenome.pl $i mm10 motifs/${i%.txt} -nomotif -bg bg/${i%_Blue*}.bed ; done
for i in grouped_peaks/*Red* ; do findMotifsGenome.pl $i mm10 motifs/${i%.txt} -nomotif -bg bg/${i%_Red*}.bed ; done
for i in grouped_peaks/*Green* ; do findMotifsGenome.pl $i mm10 motifs/${i%.txt} -nomotif -bg bg/${i%_Green*}.bed ; done
for i in grouped_peaks/*Purple* ; do findMotifsGenome.pl $i mm10 motifs/${i%.txt} -nomotif -bg bg/${i%_Purple*}.bed ; done
