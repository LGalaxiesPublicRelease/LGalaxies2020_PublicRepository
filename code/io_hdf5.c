#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"
#include "io_hdf5.h"

#ifdef HDF5_OUTPUT

void setup_hdf5_field_types(void) {

int field_ndim;
int dimProd;

// Make sure to free fields and their handles before setting-up a new HDF5 file for the next treefile:
if (field_types != NULL) {
	cleanup_hdf5_fields();
}

rowsize=0;

output_offsets=(size_t*)calloc(nfields,sizeof(output_offsets));
output_sizes=(size_t*)calloc(nfields,sizeof(output_sizes));
field_types=(hid_t*)calloc(nfields,sizeof(field_types));
is_array_type = (int*)calloc(nfields, sizeof(int));

printf("\n nfields=%d\n",nfields);
 
//This sets the types and offsets
for(ifield=0;ifield<nfields;ifield++){
 
    field_ndim=0;
    dimProd=1;
    output_offsets[ifield]=rowsize;

    // Find out how many dimensions our arrays have
    // Note: the ordering of the arrays dimensions is important
#ifdef H2_AND_RINGS
    if (flagRings[ifield]>0) {
	dims[field_ndim]=RNUM;
	dimProd*=RNUM;
	field_ndim++;
    }
#endif
#ifdef STAR_FORMATION_HISTORY	
    if (flagSFH[ifield]>0) {
	dims[field_ndim]=SFH_NBIN;
	dimProd*=SFH_NBIN;
	field_ndim++;
    }
#endif
#ifdef DETAILED_METALS_AND_MASS_RETURN
    if (flagMetals[ifield]>0) {
    	dims[field_ndim]=NUM_METAL_CHANNELS;
    	dimProd*=NUM_METAL_CHANNELS;
    	field_ndim++;
        }
    if (flagElements[ifield]>0) {
	dims[field_ndim]=NUM_ELEMENTS;
	dimProd*=NUM_ELEMENTS;
	field_ndim++;
    }
#endif
#ifdef DETAILED_DUST
    if (flagDustColdGasRates[ifield]>0) {
	dims[field_ndim]=NUM_COLDGAS_DUST_RATES;
	dimProd*=NUM_COLDGAS_DUST_RATES;
	field_ndim++;
    }
    if (flagDustHotGasRates[ifield]>0) {
	dims[field_ndim]=NUM_HOTGAS_DUST_RATES;
	dimProd*=NUM_HOTGAS_DUST_RATES;
	field_ndim++;
    }
#endif
    if (flag3[ifield]>0) {
	dims[field_ndim]=3;
	dimProd*=3;
	field_ndim++;
    }
#ifndef LITE_OUTPUT
    if (flagMag[ifield]>0) {
	dims[field_ndim]=NMAG;
	dimProd*=NMAG;
	field_ndim++;
    }	
#endif

    // Determine the basic types and their sizes
    // Structs are implemented by increasing the array dimensions by 1
    if(types[ifield]=='i'){
	field_type=H5T_NATIVE_INT;
	output_size=sizeof(int);
    }
    else if(types[ifield]=='f'){
	field_type=H5T_NATIVE_FLOAT;
	output_size=sizeof(float);
    }
    else if(types[ifield]=='l'){
	field_type=H5T_NATIVE_LLONG;
	output_size=sizeof(long long);
    }
    else terminate("Unknown field type\n");
    
#ifdef DEBUG_HDF5
	int idim; 
	printf("%s[dims]=",field_names[ifield]);
	for(idim=0;idim<field_ndim;idim++) printf("%d,",(int)dims[idim]);
	printf("\n");
	printf("output_size = %d\n",(int)output_size);
#endif

	//is_array_type = (int*)calloc(nfields, sizeof(int));

    if (field_ndim==0) {
	field_types[ifield]=field_type;
	output_sizes[ifield]=output_size;
	is_array_type[ifield]=0; //record that this is a default HDF5 type
    }
    else {
	field_types[ifield]=H5Tarray_create(field_type,field_ndim,dims);
	output_sizes[ifield]=dimProd*output_size;
	is_array_type[ifield]=1; //record that this is a custom HDF5 type
    }
    rowsize+=output_sizes[ifield];

}
}


void open_hdf5_file(int filenr, int n) {
    char output_file[80];
#ifdef GALAXYTREE
    sprintf(output_file, "%s/%s_galtree_%d.h5", OutputDir, FileNameGalaxies, filenr);
#else
    sprintf(output_file,"%s/%s_z%1.2f_%d.h5", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr); //include redshift in output file name
#endif
    file_id[n] = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}


void create_hdf5_table(int n){

  //printf("\nEntering create_hdf5_table");

  char table_name[20];
  const hsize_t nrecords=0 ; // Don't write any records on table creation

  output_size=sizeof(struct GALAXY_OUTPUT);

  struct GALAXY_OUTPUT galaxy_output;

#ifdef GALAXYTREE
  sprintf(table_name,"GalTree");
  printf("\nio_hdf5.c: Making galtree table\n");
#else
  sprintf(table_name,"%d",ListOutputSnaps[n]);
  printf("\nio_hdf5.c: Making table for snapshot %s\n",table_name);
#endif
#ifdef DEBUG_HDF5
  printf("table_name=%s\n",table_name);
  printf("file_id=%d\n",(int)file_id[n]);
  printf("nfields=%d\n",nfields);
  printf("output_size=%d\n",(int)output_size);
  int ifield;
  for (ifield=0;ifield<nfields;ifield++) {
      printf("field_names[%d]=%s\n",ifield,field_names[ifield]);
      printf("output_offsets[%d]=%d\n",ifield,(int)output_offsets[ifield]);
      printf("field_types[%d]=%d\n",ifield,(int)field_types[ifield]);
  }
#endif

  H5TBmake_table(table_name, file_id[n], table_name, nfields, nrecords,
           output_size, field_names, output_offsets, field_types,
           chunk_size, fill_data, COMPRESS, &galaxy_output);

  //printf("Leaving create_hdf5_table\n");
}


void hdf5_append_data(int n, struct GALAXY_OUTPUT * galaxy_output,int nrecords_app){

  char table_name[20];
  // Create the tablename from the snapshots
#ifdef GALAXYTREE
  sprintf(table_name,"GalTree");
#else
  sprintf(table_name,"%d",ListOutputSnaps[n]);
#endif

  H5TBappend_records(file_id[n],table_name,nrecords_app, output_size, output_offsets, output_sizes,galaxy_output);

}


void hdf5_close(int n){
  //Write tables from the input files, use: write_input_table(directory,filename)
#ifdef DEBUG_HDF5
  printf("Closing HDF5_file...\n");
#endif
  write_input_table("My_Makefile_options", n);
  write_input_table(inputFile, n);
#ifdef DEBUG_HDF5
  printf("wrote_input_table inputFile\n");
#endif
  write_prop_table(n);
#ifdef DEBUG_HDF5
  printf("wrote_prop_table\n");
#endif

  printf("\nClosing hdf5 file \n");
  H5Fclose(file_id[n]);
#ifdef DEBUG_HDF5
  printf("HDF5 file closed\n");
#endif
}


void cleanup_hdf5_fields(void){
  int ifield;
  for (ifield = 0; ifield < nfields; ifield++) {
	  if (is_array_type[ifield]) {
		  H5Tclose(field_types[ifield]);   // release custom array-type handles too
	  }
  }
  free(output_offsets);
  free(field_types);
  free(output_sizes);
  free(is_array_type);
  output_offsets = NULL;
  field_types    = NULL;
  output_sizes   = NULL;
  is_array_type = NULL;
}


void write_prop_table(int n){

  //Setup a Prop table to give a brief description about the data
  struct PropTable{
    char Name[HDF5_STRING_SIZE];
    char Units[HDF5_STRING_SIZE];
    char Description[HDF5_STRING_SIZE];
  }proptable;

  size_t proptable_output_size = sizeof(struct PropTable);

  hsize_t PropTable_nfields=3;
  const char* proptable_field_names[3]={"Field","Units","Description"};

  // Sets a string type containing CSIZE characters
  hid_t string_type;
  size_t proptable_output_offsets[3]={HOFFSET(struct PropTable,Name),
				      HOFFSET(struct PropTable,Units),
				      HOFFSET(struct PropTable,Description)};
  hid_t   proptable_field_types[3];
  size_t proptable_output_sizes[3]={sizeof(proptable.Name),
				    sizeof(proptable.Units),
				    sizeof(proptable.Description)};

  struct PropTable prop_table[127];

  string_type=H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type,HDF5_STRING_SIZE);
  proptable_field_types[0]=string_type;
  proptable_field_types[1]=string_type;
  proptable_field_types[2]=string_type;

  H5TBmake_table("Prop Table", file_id[n], "Prop Table",PropTable_nfields,NRECORDS,
		 proptable_output_size ,proptable_field_names, proptable_output_offsets, proptable_field_types,
		 chunk_size, fill_data, COMPRESS, &prop_table  );

  FILE *fp; // The file pointer 
  char line[HDF5_STRING_SIZE];
  char *token;
  const char d[2]=",";
  char filename[30]="input/hdf5_field_props.txt";

    fp=fopen(filename,"r");

  if(fp==NULL){
    printf("could not open the file %s \n", filename);
    exit(1);
  }
  // Put data into the struct and append to the table

#ifdef DEBUG_HDF5
  int icount=0;
#endif
  while(fgets(line,sizeof(line),fp)){
    token=strtok(line,d);
    if (token!=NULL) 
	strcpy(prop_table->Name,token);
    else
	strcpy(prop_table->Name," ");
    token=strtok(NULL,d);
    if (token!=NULL) 
	strcpy(prop_table->Units,token);
    else
	strcpy(prop_table->Units," ");
    token=strtok(NULL,d);
    if (token!=NULL) 
	strcpy(prop_table->Description,token);
    else
	strcpy(prop_table->Description," ");
#ifdef DEBUG_HDF5
    printf("CP3.%d\n finished prop_table assignment",icount);
    printf("prop_table->Name = %s\n",prop_table->Name);
    printf("prop_table->Description = %s\n",prop_table->Description);
    printf("prop_table->Units = %s\n",prop_table->Units);
#endif
    H5TBappend_records(file_id[n],"Prop Table",1, proptable_output_size, proptable_output_offsets, proptable_output_sizes,&prop_table);
#ifdef DEBUG_HDF5
    icount++;
#endif
  } 

  H5Tclose( string_type ); 

}


void write_input_table(char *filename, int n){

  //Setup a Prop table to give a brief description about the data
  struct InputTable{
    char Input[2048];
  }inputtable;

  char table_name[HDF5_STRING_SIZE];

  // Table filename after last slash (otherwise HDF5 crashes)
  char* delim;
  strcpy(table_name,filename);
  delim = strchr(filename,'/');
  while (delim!=NULL) {
      sprintf(table_name,"%s",++delim);
      delim = strchr(table_name,'/');
  }

  size_t inputtable_output_size = sizeof(struct InputTable);

  hsize_t InputTable_nfields=1;
  const char* inputtable_field_names[1]={filename};

  // Sets a string typ econtaining CSIZE characters
  hid_t string_type;
  size_t inputtable_output_offsets[1]={HOFFSET(struct InputTable,Input)};
  hid_t   inputtable_field_types[1];
  size_t inputtable_output_sizes[1]={sizeof(inputtable.Input)};

  struct InputTable input_table[1];

  // struct InputTable *inputtable;
  string_type=H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type,2048);
  inputtable_field_types[0]=string_type;

  H5TBmake_table(table_name, file_id[n], table_name,InputTable_nfields,NRECORDS,inputtable_output_size ,inputtable_field_names, inputtable_output_offsets, inputtable_field_types,chunk_size, fill_data, COMPRESS, &input_table  );

  FILE *fp; // The file pointer 
  char line[2048];

  fp=fopen(filename,"r");
  if(fp==NULL){
    printf("could not open the file %s \n", filename);
    exit(1);
  }
  // Put data into the struct and append to the table

 
  while(fgets(line,sizeof(line),fp)){
    strcpy(input_table->Input,line);
    H5TBappend_records(file_id[n],table_name,1, inputtable_output_size, inputtable_output_offsets, inputtable_output_sizes,&input_table);
  } 

  H5Tclose( string_type ); 
   

}


#endif //HDF5_OUTPUT
