function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s) { return rtrim(ltrim(s)); }

function getid(don) {
	chaind=substr(don,1,1);
	seriald=substr(don,2,4)+0;
	icoded=substr(don,6,1);
	idd = chaind " " seriald " " icoded;
	return idd;
}

function get_hb_no(no) {
	no_s = no;
	for(i=1;i<=6-length(no);i++)
		no_s = " " no_s;
	return no_s;
}

function set_dict(dict, t) {
	for(i in t)
		dict[t[i]]++;
}

BEGIN {
	dna_str = "DT,DA,DC,DG";
	split(dna_str,x,",");
	set_dict(dna_dict,x);
	
	rna_str = "A,C,G,U";
	#rna_str = "A,C,G,U,A23,A2L,A2M,A39,A3P,A44,A5O,A6A,A7E,A9Z,ADI,ADP,AET,AMD,AMO,AP7,AVC,MA6,MAD,MGQ,MIA,MTU,M7A,26A,2MA,6IA,6MA,6MC,6MP,6MT,6MZ,6NW,F3N,N79,RIA,V3L,ZAD,31H,31M,7AT,O2Z,SRA,00A,45A,8AN,LCA,P5P,PPU,PR5,PU,T6A,TBN,TXD,TXP,12A,1MA,5FA,A6G,E6G,E7G,EQ4,IG,IMP,M2G,MGT,MGV,MHG,QUO,YG,YYG,23G,2EG,2MG,2SG,B8K,B8W,B9B,BGH,N6G,RFJ,ZGU,7MG,CG1,G1G,G25,G2L,G46,G48,G7M,GAO,GDO,GDP,GH3,GNG,GOM,GRB,GTP,KAG,KAK,O2G,OMG,8AA,8OS,LG,PGP,P7G,TPG,TG,XTS,102,18M,1MG,A5M,A6C,E3C,IC,M4C,M5M,6OO,B8Q,B8T,B9H,JMH,N5M,RPC,RSP,RSQ,ZBC,ZCY,73W,C25,C2L,C31,C43,C5L,CBV,CCC,CH,CSF,OMC,S4C,4OC,LC,LHH,LV2,PMT,TC,10C,1SC,5HM,5IC,5MC,A6U,IU,I4U,MEP,MNU,U25,U2L,U2P,U31,U34,U36,U37,U8U,UAR,UBB,UBD,UD5,UPV,UR3,URD,US5,UZR,UMO,U23,2AU,2MU,2OM,B8H,FHU,FNU,F2T,RUS,ZBU,3AU,3ME,3MU,3TD,70U,75B,CNU,OMU,ONE,S4U,SSU,SUR,4SU,85Y,DHU,H2U,LHU,PSU,PYO,P4U,T31,125,126,127,1RN,5BU,5FU,5MU,9QV,5GP";
	split(rna_str,y,",");
	set_dict(rna_dict,y);
	
	protein_str = "ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,PRO,SER,THR,TRP,TYR,VAL";
	split(protein_str,z,",");
	set_dict(protein_dict,z);
}

length($0)==75 && $NF !~ /[0-9]+/

length($0)==75 && $NF ~ /[0-9]+/ {
	donor=substr($0,1,9);
	acceptor=substr($0,15,9);
	
	chain_d = trim(substr(donor,1,1));
	if (chain_d=="-")
		chain_d = "A";
	chain_a = trim(substr(acceptor,1,1));
	if (chain_a=="-")
		chain_a = "A";
	
	name_d = trim(substr(donor,7,3));
	name_a = trim(substr(acceptor,7,3));
	
	if ((name_d != "HOH" && name_a != "HOH") && (length(rna_chains)==0 || (length(rna_chains)>0 && (index(rna_chains,chain_d)>0 || index(rna_chains,chain_a)>0))) && 
		(rna_dict[name_d]>0 || rna_dict[name_a]>0) && !((protein_dict[name_d]>0 && protein_dict[name_a]>0) || (dna_dict[name_d]>0 && dna_dict[name_a]>0) ||
		 (rna_dict[name_d]>0 && rna_dict[name_a]>0))) {
		no=++count;
		print substr($0,1,length($0)-6) get_hb_no(no);
	}
}