#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

///////////////////////// Define paths //////////////////////////////////

char *link_vasp = "srun /opt/vasp/vasp.5.4.1_new/bin/vasp_std >stdout";
char *link_pot  = "/opt/vasp/Potential/potpaw_PBE.54";

/////////////////////Operation Section///////////////////////////////////

void write_incar_opt(char filename[300], double sigma, double magmom[8]) {
    FILE *fo;
    fo = fopen(filename, "wt");
    fprintf(fo, "# SCF Job\n");
	fprintf(fo, "#===========================================\n\n");
	
	fprintf(fo, "# General Setup\n");
	fprintf(fo, "  System = OPT Job       	# Calculation Title\n");
	fprintf(fo, "  ENCUT  = 500       		# Kinetic Energy Cutoff in eV\n");
	fprintf(fo, "  ISTART = 0         	    # Job: 0-new  1-cont  2-samecut\n");
	fprintf(fo, "  ICHARG = 2               # initial charge density: 1-file 2-atom 10-cons 11-DOS\n");
	fprintf(fo, "  ISPIN  = 2\n");
	fprintf(fo, "  ISIF   = 3\n");
	fprintf(fo, "  PREC   = Accuracy\n");
	fprintf(fo, "  IBRION = 2\n");
	fprintf(fo, "  POTIM  = 0.1\n");
	fprintf(fo, "  MAGMOM = 1*%lf 1*%lf 1*%lf 1*%lf\n", magmom[0], magmom[1], magmom[2], magmom[3]);
	fprintf(fo, "  LORBIT = 11\n");
	fprintf(fo, "  LWAVE   = .FALSE.\n");
	fprintf(fo, "  ADDGRID = .TRUE.\n\n");
	
	fprintf(fo, "# Electronic Relaxation (SCF)\n");
	fprintf(fo, "  NSW    =     300\n");
    fprintf(fo, "  NELMDL =     -6      # Number of delayed ELM steps\n\n");
	fprintf(fo, "  NELM   =     300     # Number of ELM steps\n");
	fprintf(fo, "  EDIFF  =     1.0E-7  # Stopping criteria for ESC\n");
	fprintf(fo, "  EDIFFG =     1.0E-6\n");
	fprintf(fo, "  ALGO   =     Fast    # Electronic algorithm minimization\n");
	fprintf(fo, "  SIGMA  =     %lf     # Broadening in eV\n", sigma);
	fprintf(fo, "  ISMEAR = 0           # Partial Occupancies for Each Orbital\n");
	
	fprintf(fo, "# Parallelization Scheme\n");
    fprintf(fo, "  LPLANE = .TRUE.\n");
    fprintf(fo, "  NPAR   = 4\n\n");	
	
    fclose(fo);
}

void write_incar_scf(char filename[300], double sigma, double magmom[8]) {
    FILE *fo;
    fo = fopen(filename, "wt");
    fprintf(fo, "# SCF Job\n");
	fprintf(fo, "#===========================================\n\n");
	
	fprintf(fo, "# General Setup\n");
	fprintf(fo, "  System = OPT Job       	# Calculation Title\n");
	fprintf(fo, "  ENCUT  = 500       		# Kinetic Energy Cutoff in eV\n");
	fprintf(fo, "  ISTART = 0         	    # Job: 0-new  1-cont  2-samecut\n");
	fprintf(fo, "  ICHARG = 2               # initial charge density: 1-file 2-atom 10-cons 11-DOS\n");
	fprintf(fo, "  ISPIN  = 2\n");
	fprintf(fo, "  ISIF   = 3\n");
	fprintf(fo, "  PREC   = Accuracy\n");
	fprintf(fo, "  IBRION = -1\n");
	fprintf(fo, "  POTIM  = 0.1\n");
	fprintf(fo, "  MAGMOM = 1*%lf 1*%lf 1*%lf 1*%lf\n", magmom[0], magmom[1], magmom[2], magmom[3]);
	fprintf(fo, "  LORBIT = 11\n");
	fprintf(fo, "  LWAVE   = .FALSE.\n");
	fprintf(fo, "  ADDGRID = .TRUE.\n\n");
	
	fprintf(fo, "# Electronic Relaxation (SCF)\n");
	fprintf(fo, "  NSW    =     0\n");
    fprintf(fo, "  NELMDL =     -6      # Number of delayed ELM steps\n\n");
	fprintf(fo, "  NELM   =     300     # Number of ELM steps\n");
	fprintf(fo, "  EDIFF  =     1.0E-7  # Stopping criteria for ESC\n");
	fprintf(fo, "  EDIFFG =     1.0E-6\n");
	fprintf(fo, "  ALGO   =     Fast    # Electronic algorithm minimization\n");
	fprintf(fo, "  SIGMA  =     %lf     # Broadening in eV\n", sigma);
	fprintf(fo, "  ISMEAR = 0           # Partial Occupancies for Each Orbital\n");
	fprintf(fo, "  NEDOS  = 3000\n");
	
	fprintf(fo, "# Parallelization Scheme\n");
    fprintf(fo, "  LPLANE = .TRUE.\n");
    fprintf(fo, "  NPAR   = 4\n\n");	
	
    fclose(fo);
}

void write_incar_opt2(char filename[300], double sigma, double magmom[8]) {
    FILE *fo;
    fo = fopen(filename, "wt");
    fprintf(fo, "# SCF Job\n");
	fprintf(fo, "#===========================================\n\n");
	
	fprintf(fo, "# General Setup\n");
	fprintf(fo, "  System = OPT Job       	# Calculation Title\n");
	fprintf(fo, "  ENCUT  = 500       		# Kinetic Energy Cutoff in eV\n");
	fprintf(fo, "  ISTART = 0         	    # Job: 0-new  1-cont  2-samecut\n");
	fprintf(fo, "  ICHARG = 2               # initial charge density: 1-file 2-atom 10-cons 11-DOS\n");
	fprintf(fo, "  ISPIN  = 2\n");
	fprintf(fo, "  ISIF   = 3\n");
	fprintf(fo, "  PREC   = Accuracy\n");
	fprintf(fo, "  IBRION = 2\n");
	fprintf(fo, "  POTIM  = 0.1\n");
	fprintf(fo, "  MAGMOM = 1*%lf 1*%lf 1*%lf 1*%lf\n", magmom[4], magmom[5], magmom[6], magmom[7]);
	fprintf(fo, "  LORBIT = 11\n");
	fprintf(fo, "  LWAVE   = .FALSE.\n");
	fprintf(fo, "  ADDGRID = .TRUE.\n\n");
	
	fprintf(fo, "# Electronic Relaxation (SCF)\n");
	fprintf(fo, "  NSW    =     300\n");
    fprintf(fo, "  NELMDL =     -6      # Number of delayed ELM steps\n\n");
	fprintf(fo, "  NELM   =     300     # Number of ELM steps\n");
	fprintf(fo, "  EDIFF  =     1.0E-7  # Stopping criteria for ESC\n");
	fprintf(fo, "  EDIFFG =     1.0E-6\n");
	fprintf(fo, "  ALGO   =     Fast    # Electronic algorithm minimization\n");
	fprintf(fo, "  SIGMA  =     %lf     # Broadening in eV\n", sigma);
	fprintf(fo, "  ISMEAR = 0           # Partial Occupancies for Each Orbital\n");
	
	fprintf(fo, "# Parallelization Scheme\n");
    fprintf(fo, "  LPLANE = .TRUE.\n");
    fprintf(fo, "  NPAR   = 4\n\n");	
	
    fclose(fo);
}

void write_incar_scf2(char filename[300], double sigma, double magmom[8]) {
    FILE *fo;
    fo = fopen(filename, "wt");
    fprintf(fo, "# SCF Job\n");
	fprintf(fo, "#===========================================\n\n");
	
	fprintf(fo, "# General Setup\n");
	fprintf(fo, "  System = OPT Job       	# Calculation Title\n");
	fprintf(fo, "  ENCUT  = 500       		# Kinetic Energy Cutoff in eV\n");
	fprintf(fo, "  ISTART = 0         	    # Job: 0-new  1-cont  2-samecut\n");
	fprintf(fo, "  ICHARG = 2               # initial charge density: 1-file 2-atom 10-cons 11-DOS\n");
	fprintf(fo, "  ISPIN  = 2\n");
	fprintf(fo, "  ISIF   = 3\n");
	fprintf(fo, "  PREC   = Accuracy\n");
	fprintf(fo, "  IBRION = -1\n");
	fprintf(fo, "  POTIM  = 0.1\n");
	fprintf(fo, "  MAGMOM = 1*%lf 1*%lf 1*%lf 1*%lf\n", magmom[4], magmom[5], magmom[6], magmom[7]);
	fprintf(fo, "  LORBIT = 11\n");
	fprintf(fo, "  LWAVE   = .FALSE.\n");
	fprintf(fo, "  ADDGRID = .TRUE.\n\n");
	
	fprintf(fo, "# Electronic Relaxation (SCF)\n");
	fprintf(fo, "  NSW    =     0\n");
    fprintf(fo, "  NELMDL =     -6      # Number of delayed ELM steps\n\n");
	fprintf(fo, "  NELM   =     300     # Number of ELM steps\n");
	fprintf(fo, "  EDIFF  =     1.0E-7  # Stopping criteria for ESC\n");
	fprintf(fo, "  EDIFFG =     1.0E-6\n");
	fprintf(fo, "  ALGO   =     Fast    # Electronic algorithm minimization\n");
	fprintf(fo, "  SIGMA  =     %lf     # Broadening in eV\n", sigma);
	fprintf(fo, "  ISMEAR = 0           # Partial Occupancies for Each Orbital\n");
	fprintf(fo, "  NEDOS  = 3000\n");
	
	fprintf(fo, "# Parallelization Scheme\n");
    fprintf(fo, "  LPLANE = .TRUE.\n");
    fprintf(fo, "  NPAR   = 4\n\n");	
	
    fclose(fo);
}

void write_kpoints(char filename[300],int k[3]) {
    FILE *fo;
    fo = fopen(filename, "wt");
    fprintf(fo, "Automatic Mesh # Generates Automatically the K-mesh\n");
    fprintf(fo, "  0\n");
    fprintf(fo, "Gamma # Automatic Distribution of the mesh\n");
	fprintf(fo, "  %d  %d  %d\n",k[0],k[1],k[2]);
    fclose(fo);
}


void write_poscar_l21(char filename[300], char X[10], char Y[10], char Z[10],  double lat) {
    FILE *fo;
    fo = fopen(filename, "a");
    fprintf(fo, "%s %s %s\n", X, Y, Z);
	fprintf(fo, "   %lf\n", lat);
	fprintf(fo, "     0.0000000000000000    0.5000000000000000    0.5000000000000000\n");
	fprintf(fo, "     0.5000000000000000    0.0000000000000000    0.5000000000000000\n");
	fprintf(fo, "     0.5000000000000000    0.5000000000000000    0.0000000000000000\n");
	fprintf(fo, "   %s   %s    %s\n", X, Y, Z);
	fprintf(fo, "    2    1     1\n");
	fprintf(fo, "Direct\n");
	fprintf(fo, "  0.7500000000000000  0.7500000000000000  0.7500000000000000\n");
	fprintf(fo, "  0.2500000000000000  0.2500000000000000  0.2500000000000000\n");
	fprintf(fo, "  0.0000000000000000  0.0000000000000000  0.0000000000000000\n");
	fprintf(fo, "  0.5000000000000000  0.5000000000000000  0.5000000000000000\n");
    fclose(fo);
}

void write_poscar_xa(char filename[300], char X[10], char Y[10], char Z[10],  double lat) {
    FILE *fo;
    fo = fopen(filename, "a");
    fprintf(fo, "%s %s %s\n", X, Y, Z);
	fprintf(fo, "   %lf\n", lat);
	fprintf(fo, "     0.0000000000000000    0.5000000000000000    0.5000000000000000\n");
	fprintf(fo, "     0.5000000000000000    0.0000000000000000    0.5000000000000000\n");
	fprintf(fo, "     0.5000000000000000    0.5000000000000000    0.0000000000000000\n");
	fprintf(fo, "   %s   %s    %s\n", X, Y, Z);
	fprintf(fo, "    2    1     1\n");
	fprintf(fo, "Direct\n");
	fprintf(fo, "  0.0000000000000000  0.0000000000000000  0.0000000000000000\n");
	fprintf(fo, "  0.2500000000000000  0.2500000000000000  0.2500000000000000\n");
	fprintf(fo, "  0.7500000000000000  0.7500000000000000  0.7500000000000000\n");
	fprintf(fo, "  0.5000000000000000  0.5000000000000000  0.5000000000000000\n");
    fclose(fo);
}

int main() {
	// initial steps for getting parameters
	int i, j;
	char filename[300];
	char pot1[300], pot2[300], pot3[300];
	char line[300];
	int num_struct;
	
	int kpoint[3][3];
	kpoint[0][0] = 3;kpoint[0][1]  = 3; kpoint[0][2]  = 3;
	kpoint[1][0] = 9;kpoint[1][1]  = 9; kpoint[1][2]  = 9;
	kpoint[2][0] = 15;kpoint[2][1] = 15; kpoint[2][2] = 15;

	double sigma[4]; sigma[0] = 0.2; sigma[1] = 0.1; sigma[2] = 0.1; sigma[3] = 0.05;
	
	// read INPUT
	FILE *fi;
    fi = fopen("INPUT", "rt");
	fgets(line,300,fi);
	sscanf(line,"%d",&num_struct);
	char elements[num_struct][3][300];
	double lat[num_struct];
	
	char type_pot[num_struct][3][300];
	
	double magmom[num_struct][8];
	
	for (i = 0;i < num_struct;i++){
		fgets(line,300,fi);
		sscanf(line,"%s %s %s %lf %s %s %s %lf %lf %lf %lf %lf %lf %lf %lf",&elements[i][0], &elements[i][1], &elements[i][2], &lat[i], &type_pot[i][0], &type_pot[i][1], &type_pot[i][2], &magmom[i][0], &magmom[i][1], &magmom[i][2], &magmom[i][3], &magmom[i][4], &magmom[i][5], &magmom[i][6], &magmom[i][7]);
	}
	
	fclose(fi);
	
	// Start calculations
	for (i = 0;i < num_struct;i++){
		
		sprintf(filename, "mkdir %s%s%s",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		//////// Calculate l21 structure //////
		sprintf(filename, "mkdir %s%s%s/l21",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		// opt
		sprintf(filename, "mkdir %s%s%s/l21/opt",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		sprintf(filename, "mkdir %s%s%s/l21/opt/data",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		// make input files
		if (strcmp(type_pot[i][0],"pv") == 0){
			sprintf(pot1,"%s/%s_pv/POTCAR",link_pot, elements[i][0]);
		}
		else if (strcmp(type_pot[i][0],"sv") == 0){
			sprintf(pot1,"%s/%s_sv/POTCAR",link_pot, elements[i][0]);
		}
		else{
			sprintf(pot1,"%s/%s/POTCAR",link_pot, elements[i][0]);
		}
		
		if (strcmp(type_pot[i][1],"pv") == 0){
			sprintf(pot2,"%s/%s_pv/POTCAR",link_pot, elements[i][1]);
		}
		else if (strcmp(type_pot[i][1],"sv") == 0){
			sprintf(pot2,"%s/%s_sv/POTCAR",link_pot, elements[i][1]);
		}
		else{
			sprintf(pot2,"%s/%s/POTCAR",link_pot, elements[i][1]);
		}
		
		if (strcmp(type_pot[i][2],"pv") == 0){
			sprintf(pot3,"%s/%s_pv/POTCAR",link_pot, elements[i][2]);
		}
		else if (strcmp(type_pot[i][2],"sv") == 0){
			sprintf(pot3,"%s/%s_sv/POTCAR",link_pot, elements[i][2]);
		}
		else{
			sprintf(pot3,"%s/%s/POTCAR",link_pot, elements[i][2]);
		}
		
		sprintf(filename, "cat %s %s %s > %s%s%s/l21/opt/POTCAR", pot1, pot2, pot3, elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		sprintf(filename, "%s%s%s/l21/opt/POSCAR",elements[i][0], elements[i][1], elements[i][2]);
		write_poscar_l21(filename,elements[i][0], elements[i][1], elements[i][2], lat[i]);
		for (j = 0;j < 3;j++){
			sprintf(filename, "%s%s%s/l21/opt/INCAR",elements[i][0], elements[i][1], elements[i][2]);
			write_incar_opt(filename,sigma[j], magmom[i]);
			sprintf(filename, "%s%s%s/l21/opt/KPOINTS",elements[i][0], elements[i][1], elements[i][2]);
			write_kpoints(filename,kpoint[j]);
			// run vasp
			sprintf(filename, "cd %s%s%s/l21/opt && %s && cd .. && cd .. && cd ..", elements[i][0], elements[i][1], elements[i][2], link_vasp);
	        system(filename);
			
			// backup data
			sprintf(filename, " cp %s%s%s/l21/opt/CONTCAR %s%s%s/l21/opt/data/POSCAR_%d",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2], kpoint[j][0]);
			system(filename);
			sprintf(filename, " cp %s%s%s/l21/opt/OUTCAR %s%s%s/l21/opt/data/OUTCAR_%d",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2], kpoint[j][0]);
			system(filename);
			sprintf(filename, " cp %s%s%s/l21/opt/OSZICAR %s%s%s/l21/opt/data/OSZICAR_%d",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2], kpoint[j][0]);
			system(filename);
			sprintf(filename, " cp %s%s%s/l21/opt/CONTCAR %s%s%s/l21/opt/POSCAR",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2]);
			system(filename);
		}
		
		// scf
		sprintf(filename, "mkdir %s%s%s/l21/scf",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		// copy files
		sprintf(filename, "cp %s%s%s/l21/opt/POSCAR %s%s%s/l21/scf/POSCAR",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		sprintf(filename, "cp %s%s%s/l21/opt/POTCAR %s%s%s/l21/scf/POTCAR",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		// write inputs
		sprintf(filename, "%s%s%s/l21/scf/INCAR",elements[i][0], elements[i][1], elements[i][2]);
		write_incar_scf(filename,sigma[3], magmom[i]);
		sprintf(filename, "%s%s%s/l21/scf/KPOINTS",elements[i][0], elements[i][1], elements[i][2]);
		write_kpoints(filename,kpoint[2]);
		// run vasp
		sprintf(filename, "cd %s%s%s/l21/scf && %s && cd .. && cd .. && cd ..", elements[i][0], elements[i][1], elements[i][2], link_vasp);
	    system(filename);
		
		//////// Calculate xa structure //////
		sprintf(filename, "mkdir %s%s%s/xa",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		// opt
		sprintf(filename, "mkdir %s%s%s/xa/opt",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		sprintf(filename, "mkdir %s%s%s/xa/opt/data",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		// make input files
		sprintf(filename, "cat %s %s %s > %s%s%s/xa/opt/POTCAR", pot1, pot2, pot3, elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		sprintf(filename, "%s%s%s/xa/opt/POSCAR",elements[i][0], elements[i][1], elements[i][2]);
		write_poscar_xa(filename,elements[i][0], elements[i][1], elements[i][2], lat[i]);
		for (j = 0;j < 3;j++){
			sprintf(filename, "%s%s%s/xa/opt/INCAR",elements[i][0], elements[i][1], elements[i][2]);
			write_incar_opt2(filename,sigma[j], magmom[i]);
			sprintf(filename, "%s%s%s/xa/opt/KPOINTS",elements[i][0], elements[i][1], elements[i][2]);
			write_kpoints(filename,kpoint[j]);
			// run vasp
			sprintf(filename, "cd %s%s%s/xa/opt && %s && cd .. && cd .. && cd ..", elements[i][0], elements[i][1], elements[i][2], link_vasp);
	        system(filename);
			
			// backup data
			sprintf(filename, " cp %s%s%s/xa/opt/CONTCAR %s%s%s/xa/opt/data/POSCAR_%d",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2], kpoint[j][0]);
			system(filename);
			sprintf(filename, " cp %s%s%s/xa/opt/OUTCAR %s%s%s/xa/opt/data/OUTCAR_%d",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2], kpoint[j][0]);
			system(filename);
			sprintf(filename, " cp %s%s%s/xa/opt/OSZICAR %s%s%s/xa/opt/data/OSZICAR_%d",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2], kpoint[j][0]);
			system(filename);
			sprintf(filename, " cp %s%s%s/xa/opt/CONTCAR %s%s%s/xa/opt/POSCAR",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2]);
			system(filename);
		}
		
		// scf
		sprintf(filename, "mkdir %s%s%s/xa/scf",elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		// copy files
		sprintf(filename, "cp %s%s%s/xa/opt/POSCAR %s%s%s/xa/scf/POSCAR",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		sprintf(filename, "cp %s%s%s/xa/opt/POTCAR %s%s%s/xa/scf/POTCAR",elements[i][0], elements[i][1], elements[i][2], elements[i][0], elements[i][1], elements[i][2]);
		system(filename);
		
		// write inputs
		sprintf(filename, "%s%s%s/xa/scf/INCAR",elements[i][0], elements[i][1], elements[i][2]);
		write_incar_scf2(filename,sigma[3], magmom[i]);
		sprintf(filename, "%s%s%s/xa/scf/KPOINTS",elements[i][0], elements[i][1], elements[i][2]);
		write_kpoints(filename,kpoint[2]);
		// run vasp
		sprintf(filename, "cd %s%s%s/xa/scf && %s && cd .. && cd .. && cd ..", elements[i][0], elements[i][1], elements[i][2], link_vasp);
	    system(filename);
		
	}		
   
return 0;
}
