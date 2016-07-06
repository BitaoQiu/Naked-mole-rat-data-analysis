for i in *.finale.DEG;
do
/share/raid12/qiubitao/bin/bin/python3 refine_deg.py $i $i.S $i.filter
done
