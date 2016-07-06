for i in */;do
cd $i;
qsub -cwd -l vf=0.7G -q ngb.q -P ngb_un run.sh;
cd ..
done
