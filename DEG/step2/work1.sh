#cp ../Others/*/*compare4 ./ 
for i in *.compare4_3;
do
grep 'HGL_' $i | cut -f1-3,5 > $i.refine
done
for i in *.compare4_3.refine;
do
/share/raid12/qiubitao/bin/bin/python3 Make_format.py $i > $i.finale
done
for i in *.finale;
do
/share/raid12/qiubitao/bin/bin/python3 Get_deg.py $i > $i.DEG
done

