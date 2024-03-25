function print_selected_residue(id, a, all) {
	for(i in a) {
		if (a[i] == id) {
			print get_res_no(++no) substr(all,6);
			break;
		}
	}
}

function get_res_no(no) {
	no_s = no;
	for(i=1;i<=5-length(no);i++)
		no_s = " " no_s;
	return no_s;
}

BEGIN {split(residues,a,",");}
$1 == "#"
$3 ~ /[A-Za-z0-9]+/ && $2 ~ /[0-9]+/ {
	id = $3 $2;
	print_selected_residue(id,a,$0);
}
