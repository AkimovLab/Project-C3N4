for j in {FSSH,IDA,mSDM}; do 
   rm $j.txt
   cat out_e-h-recom_c3n4.log | grep 'Special flag' | grep "$j" | awk '{print $14 " " $16}' >> $j.txt ; 
done
