GET_HOMOLOGUES Step-by-Step Guide
Create and activate a conda environment with required tools:
``` bash
conda create -n gethomologues_env perl blast mafft hmmer -y
conda activate gethomologues_env

``` 
Install BioPerl (required by GET_HOMOLOGUES):
``` bash
cpan App::cpanminus
cpanm Bio::Perl
``` 
Install GET_HOMOLOGUES
``` bash
git clone https://github.com/eead-csic-compbio/get_homologues.git ~/get_homologues
export PATH=$HOME/get_homologues:$PATH
echo 'export PATH=$HOME/get_homologues:$PATH' >> ~/.bashrc
``` 
Check installation:
``` bash
get_homologues.pl -h
``` 
