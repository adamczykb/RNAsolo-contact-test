BEGIN {FS="\t";}
NR==1 {
	for(i=2;i<=NF;i++) {
		split($i,a,"#");
		split(a[2],b,":");
		residues[i-1]=b[2] " " b[1];
	}
}
NR>1 {
	split($1,c,":");
	count = 0;
	for(i=2;i<=NF;i++)
		if ($i == 1) {
			selected_residues[residues[i-1]]++;
			count++;
		}
	if (count>0)
		selected_residues[c[3] " " c[2]]++;
}
END {
	for(i in selected_residues)
		print i;
}