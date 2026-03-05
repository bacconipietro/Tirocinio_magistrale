# Alias
```bash
nano ~/.bashrc

#copy this line at the bottom of the file
alias unibo="ssh pietro.bacconi@studio.unibo.it@137.204.142.152"
source ~/.bashrc

#write to run
unibo
```
# Screen

Screen working
```bash
screen -S <name screen>  //create screen de novo
screen -ls //list opne screen
screen -r <name screen> //open old screen 
```

- Ctrl+A+D  //close screen 

# zip
```bash
zip -r folder_name.zip <folder>
```

# Tree directory structure 
```bash
wget https://man.archlinux.org/man/tree.1.en
[base]
conda install -c bioconda tree
tree -d Tirocinio_magistrale/
```
