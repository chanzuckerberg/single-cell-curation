in_dir='./data/remixed/'

for i in ${in_dir}/*
do
    echo "Working with $i"
    cellxgene schema validate $i
done
    
