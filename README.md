# gpcr-gprotein-3D-comparison

This page provides the scripts and data of the analysis of Gi- and Gs-GPCR 3D complex structure, both for contact as well as RMSD analysis, from Okamato et al., Nature Structure & Molecular Biology, 2021.


A local installation of the most recent versions of Psi-blast, Clustalo and Cytoscape are required to reproduce the pipeline.
We recommend launching the python scripts via anaconda environment (Python3) after having installed the most recent version of the following libraries:

-numpy

-scipy

-biopython

Steps to perform 3D superimposition and RMSD clustering via bash command line:

cd 3D_fit/

./split_fasta.py < Gprot_chains.txt > Gprot.fa

./split_fasta.py < GPCR_chains.txt > GPCR.fa

clustalo -hmm-in=7tm_1.hmm -i GPCR.fa -o GPCR_ali.fa --outfmt=fasta

clustalo -hmm-in=G-alpha.hmm -i Gprot.fa -o Gprot_ali.fa --outfmt=fasta

./msa_consensus_mappdb.py Gprot_ali.fa > Gprot_pdbs_consensus.txt 

./msa_consensus_mappdb.py GPCR_ali.fa > GPCR_pdbs_consensus.txt 

./gen_bioalign.py > run_bioalign.sh 

./run_bioalign.sh

Script to generate the plots for the figure

./rmsd_clustering.py

./rmsd_violin_plot.py


Steps to perform contact network analysis via bash command line and Cytoscape:


cd contact/

./run_get_contact.sh

./contact_analysis_comp.py ../Gi_pdbs.txt > Gi_contacts.log

./contact_analysis_comp.py ../Gs_pdbs.txt > Gs_contacts.log


Script to generate the nodes.txt and links.txt file to load and visualize the contact network through Cytoscape

./Gi_Gs_contact_net_comparison.py



