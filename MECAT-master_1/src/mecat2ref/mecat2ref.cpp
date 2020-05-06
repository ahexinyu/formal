#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/time.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include "output.h"
#include "../common/defs.h"
#define RM 100000
#define big_size 1000000
#define FM 1000000
#define MAXSTR 1000000000
#define split_le 15000

static const char* prog_name = NULL;
static const int kDefaultNumCandidates = 10;
static const int kDefaultNumOutput = 10;
static int num_candidates;
static int num_output;
static int output_format = FMT_REF;
static const int kDefaultOutputFormat = FMT_REF;
static int tech;
static const int kDefaultTech = TECH_PACBIO;




typedef struct
{
	const char* reads;
	const char* reference;
	const char* wrk_dir;
	const char* output;
    const char* refoutput;
	int         num_cores;
	int			num_candidates;
	int			num_output;
	int			output_format;
	int 		tech;
} meap_ref_options;

void init_meap_ref_options(meap_ref_options* options)
{
	options->reads = NULL;
	options->reference = NULL;
	options->wrk_dir = NULL;
	options->output = NULL;
    options->refoutput=NULL;
	options->num_cores = 1;
	options->num_candidates = kDefaultNumCandidates;
	options->num_output = kDefaultNumOutput;
	options->output_format = kDefaultOutputFormat;
	options->tech = kDefaultTech;
}

void print_usage()
{
	fprintf(stderr, "\n\n");
	fprintf(stderr, "usage:\n");
	fprintf(stderr, "%s [-d reads] [-r reference] [-o output][-p refoutput] [-w working dir] [-t threads]", prog_name);
	fprintf(stderr, "\n\n");
	fprintf(stderr, "options:\n");
	fprintf(stderr, "-d <string>\treads file name\n");
	fprintf(stderr, "-r <string>\treference file name\n");
	fprintf(stderr, "-o <string>\toutput file name\n");
    fprintf(stderr, "-p <string>\trefoutput file name\n");
	fprintf(stderr, "-w <string>\tworking folder name, will be created if not exist\n");
	fprintf(stderr, "-t <integer>\tnumber of cput threads\n\t\tdefault: 1\n");
	fprintf(stderr, "-n <integer>\tnumber of of candidates for gap extension\n\t\tdefault: %d\n", kDefaultNumCandidates);
	fprintf(stderr, "-b <integer>\toutput the best b alignments\n\t\tdefault: %d\n", kDefaultNumOutput);
	fprintf(stderr, "-m <0/1/2>\toutput format: 0 = ref, 1 = m4, 2 = sam\n\t\tdefault: %d\n", kDefaultOutputFormat);
	fprintf(stderr, "-x <0/1>\tsequencing technology: 0 = pacbio, 1 = nanopore\n\t\tdefault: %d\n", kDefaultTech);
}

int
param_read_t(int argc, char* argv[], meap_ref_options* options)
{
	int opt_char;
	char err_char;
	opterr = 0;
	int ret = 1;
	
	init_meap_ref_options(options);
	while((opt_char = getopt(argc, argv, "d:r:w:o:p:t:n:b:m:x:")) != -1)
	{
		switch(opt_char)
		{
			case 'd':
				options->reads = optarg;
				break;
			case 'r':
				options->reference = optarg;
				break;
			case 'w':
				options->wrk_dir = optarg;
				break;
			case 'o':
				options->output = optarg;
				break;
            case 'p':
                options->refoutput = optarg;
                break;
			case 't':
				options->num_cores = atoi(optarg);
				break;
			case 'n':
				options->num_candidates = atoi(optarg);
				break;
			case 'b':
				options->num_output = atoi(optarg);
				break;
			case 'm':
				options->output_format = atoi(optarg);
				break;
			case 'x':
				if (optarg[0] == '0') {
					options->tech = TECH_PACBIO;
				} else if (optarg[0] == '1') {
					options->tech = TECH_NANOPORE;
				} else {
					ERROR("Invalid argument to option 'x': %s\n", optarg);
				}
				break;
			case ':':
				err_char = (char)optopt;
				fprintf(stderr, "Error: unrecogised option \'%c\'\n", err_char);
				ret = -1;
				return ret;
			case '?':
				err_char = (char)optopt;
				fprintf(stderr, "Error: argument to option \'%c\' is missing!\n", err_char);
				ret = -1;
				return ret;	
		}
	}
	
	const char* options_err_msg = NULL;
	if (!options->reads)
		options_err_msg = "dataset must be specified";
	else if (!options->reference)
		options_err_msg = "reference must be specified";
	else if (!options->output)
		options_err_msg = "output must be specified";
    else if (!options->refoutput)
        options_err_msg = "refoutput must be specified";
	else if (!options->wrk_dir)
		options_err_msg = "working directory must be specified";
	else if (options->num_cores < 1)
		options_err_msg = "cpu cores must be > 0";
	else if (options->num_candidates < 1)
		options_err_msg = "candidates must be > 0";
	else if (options->num_output < 1)
		options_err_msg = "output alignments must be > 0";
	if (options_err_msg)
	{
		fprintf(stderr, "Error: %s\n", options_err_msg);
		ret = -1;
		return ret;
	}
	
	if (options->num_output > options->num_candidates)
	{
		fprintf(stderr, "warning: number of output (%d) is greater than number of candidates (%d), we reset it to %d",
			options->num_output, options->num_candidates, options->num_candidates);
		options->num_output = options->num_candidates;
	}
	
	DIR* dirptr = opendir(options->wrk_dir);
    if (dirptr == NULL)
    {
        int t = mkdir(options->wrk_dir, S_IRWXU);
        if (t == -1)
        {
            fprintf(stderr, "Fail to create folder %s!\n", options->wrk_dir);
            ret = -1;
        }
    }
    else
    {
        closedir(dirptr);
    }

    return ret;
}

int str2num(char *str)
{

    int sum = 0, i;
    char ch;
    i = 0;
    ch = str[i];
    while (ch != '\0')
    {
        sum = sum * 10 + (ch - '0');
        i++;
        ch = str[i];
    }
    return (sum);
}
static long get_file_size(const char *path)
{
    long filesize = -1;
    struct stat statbuff;
    if(stat(path, &statbuff) < 0)
    {
        return filesize;
    }
    else
    {
        filesize = statbuff.st_size;
    }
    return filesize;
}


static unsigned short atcttrans(char c)
{
    if(c=='A')return 0;
    else if(c=='T')return 1;
    else if(c=='C')return 2;
    else if(c=='G')return 3;
    else return 4;
}
static long sumvalue_x(int *intarry,int count)
{//countin countindex
    long i,sumval=0;
    for(i=0; i<count; i++)
    {
        
        if(intarry[i]>0&&intarry[i]<129)sumval=sumval+intarry[i];
        else if(intarry[i]>128)intarry[i]=0;
    }
    //printf("intarry is %d\n",intarry[149]);
    return(sumval);
}
int chang_fastqfile(const char *fastaq, const char *fenfolder)
{
    FILE *fp, *ot;
    int kk = 0,read_len;
    char read_name[1000], onedata[RM], buff1[1000], buff2[RM], *fq, tempstr[200], *oq,ch;
    sprintf(tempstr, "%s/0.fq", fenfolder);
    fp = fopen(fastaq, "r");
    ot = fopen(tempstr, "w");
    fq = (char *)malloc(100000000);
    setvbuf(fp, fq, _IOFBF, 100000000);
    oq = (char *)malloc(100000000);
    setvbuf(ot, oq, _IOFBF, 100000000);
	int num_read_items;
    ch=getc(fp);
    if(ch=='>')
    {
        kk=0;
        for (; ch!=EOF; ch=getc(fp))
        {
            if(ch=='>')
            {
                num_read_items = fscanf(fp,"%[^\n]s",read_name);
				assert(num_read_items = 1);
                if(kk>0)
                {
                    onedata[read_len]='\0';
                    fprintf(ot, "%d\t\%d\t%s\n", kk-1,read_len,onedata);
                    read_len=0;
                    kk++;
                }
                else
                {
                    kk++;
                    read_len=0;
                }
            }
            else if(ch!='\n'&&ch!='\r')onedata[read_len++]=ch;
        }
        onedata[read_len]='\0';
        fprintf(ot, "%d\t\%d\t%s\n", kk-1,read_len,onedata);
    }
    else
    {
        fseek(fp, 0L, SEEK_SET);
        while (fscanf(fp, "%[^\n]s", read_name) != EOF && fscanf(fp, "%s\n", onedata) != EOF && fscanf(fp, "%[^\n]s", buff1) != EOF && fscanf(fp, "%s\n", buff2) != EOF)
        {
            read_len=strlen(onedata);
            fprintf(ot, "%d\t%d\t%s\n", ++kk,read_len,onedata);
        }
    }
    fclose(fp);
    fclose(ot);
    free(fq);
    free(oq);
    return (kk);

}
//*******更改reference文件格式****************
 int change_ref_fq(const char *filepath,const char *outpath){
    
    FILE *fp;FILE *ot;int ref_len;long i;
    char tempstr[200];char *onedata;char ch;char *fq,*oq;char *buff1, *buff2;
    fp=fopen(filepath,"r");
    char refname[1000];
    sprintf(tempstr,"%s/ref.fq",outpath);
    ot=fopen(tempstr,"w");
     
     onedata=(char *)malloc(600000000);
    /*buff2=(char *)malloc(100000000);
    buff1=(char *)malloc(100000000);
    fq = (char *)malloc(1000000000);
    setvbuf(fp, fq, _IOFBF, 1000000000);
    oq = (char *)malloc(1000000000);
    setvbuf(ot, oq, _IOFBF, 1000000000);*/
    int num_read_items;
    ch=getc(fp);int kk=0;
     
   if(ch=='>')
    {
        kk=0;
        for (; ch!=EOF; ch=getc(fp))
        {
            if(ch=='>')
            {
                num_read_items = fscanf(fp,"%[^\n]s",refname);
                assert(num_read_items = 1);
                if(kk>0)
                {
                    onedata[ref_len]='\0';
                    fprintf(ot, "%d\t\%d\t%s\n", kk-1,ref_len,onedata);
                    ref_len=0;
                    kk++;
                }
                else
                {
                    kk++;
                    ref_len=0;
                }
            }
            else if(ch!='\n'&&ch!='\r')onedata[ref_len++]=ch;
        }
        onedata[ref_len]='\0';
        fprintf(ot, "%d\t\%d\t%s\n", kk-1,ref_len,onedata);
       
    }
    else
    {
        fseek(fp, 0L, SEEK_SET);
        while (fscanf(fp, "%[^\n]s", refname) != EOF && fscanf(fp, "%s\n", onedata) != EOF && fscanf(fp, "%[^\n]s", buff1) != EOF && fscanf(fp, "%s\n", buff2) != EOF)
        {
            ref_len=strlen(onedata);
            fprintf(ot, "%d\t%d\t%s\n", ++kk,ref_len,onedata);
        }
        
    }
    
    /*fclose(fp);
    fclose(ot);
    free(fq);
    free(oq);*/
    return (kk);
     //return 1;
}


int firsttask(int argc, char *argv[])
{
    struct timeval tpstart, tpend;
    float timeuse;
    FILE *fp;
	meap_ref_options* options = (meap_ref_options*)malloc(sizeof(meap_ref_options));
   
	int flag = param_read_t(argc, argv, options);
	if (flag == -1) { print_usage(); exit(1); }
  
    int readcount = chang_fastqfile(options->reads, options->wrk_dir);
    
    int ref_count = change_ref_fq(options->reference,options->wrk_dir);//*********
    
    /*char kkkkk[1024];
    sprintf(kkkkk, "config.txt");
    FILE* fileout = fopen(kkkkk, "w");
    fprintf(fileout, "%s\n%s\n%s\n%s\n%s\n%d\t%d\n%d\n", options->wrk_dir, options->reference, options->reads, options->output,options->refoutput,options->num_cores, readcount,ref_count);
    fclose(fileout);
	int corenum = options->num_cores;
	num_candidates = options->num_candidates;
	num_output = options->num_output;
	output_format = options->output_format;
	tech = options->tech;
	free(options);
    return (corenum);*/
    //return 1;
    
}

int cmp_temp_result_ptr(const void* a, const void* b)
{
	TempResult* pa = *(TempResult**)a;
	TempResult* pb = *(TempResult**)b;
	if (pa->aln_size > pb->aln_size) return -1;
	if (pa->aln_size == pb->aln_size) return 0;
	return 1;
}

int get_chr_id(fastaindexinfo* chr_idx, const int num_chr, const long offset)
{
	int left = 0, right = num_chr, mid = 0;
	while (left < right)
	{
		mid = (left + right) >> 1;
		if (offset >= chr_idx[mid].chrstart)
		{
			if (mid == num_chr - 1) break;
			if (offset < chr_idx[mid + 1].chrstart) break;
			left = mid + 1;
		}
		else
		{
			right = mid;
		}
	}
	return mid;
}

void
filter_contained_results(TempResult** pptr, const int num_results, int* valid)
{
	int i = 0, j;
	int b1, e1, id1, b2, e2, id2;
	
	for (i = 0; i < num_results; ++i) valid[i] = 1;
	for (i = 0; i < num_results - 1; ++i)
	{
		if (!valid[i]) continue;
		b1 = pptr[i]->qb;
		e1 = pptr[i]->qe;
		id1 = pptr[i]->read_id;
		for (j = i + 1; j < num_results; ++j)
		{
			if (!valid[j]) continue;
			b2 = pptr[j]->qb;
			e2 = pptr[j]->qe;
			id2 = pptr[j]->read_id;
			if (id2 != id1) continue;
			if (pptr[i]->read_dir != pptr[j]->read_dir) continue;
			if (b2 >= b1 && e2 <= e1) valid[j] = 0;
		}
	}
}

void
output_query_results(fastaindexinfo* chr_idx, const int num_chr, TempResult** pptr, const int num_results, FILE* out)
{
	//int valid[num_results];
	//qsort(pptr, num_results, sizeof(TempResult*), cmp_temp_result_ptr);
	//filter_contained_results(pptr, num_results, valid);
	const int n = num_results;
	int output_cnt = 0;
	int i;
	for (i = 0; i < n; ++i)
	{
		//if (!valid[i]) continue;
		int sid = get_chr_id(chr_idx, num_chr, pptr[i]->sb);
		output_one_result(pptr[i]->read_id,
						  chr_idx[sid].chrname,
						  pptr[i]->read_dir,
						  pptr[i]->qb,
						  pptr[i]->qe,
						  pptr[i]->qs,
						  pptr[i]->vscore,
						  pptr[i]->sb - chr_idx[sid].chrstart,
						  pptr[i]->se - chr_idx[sid].chrstart,
						  chr_idx[sid].chrsize,
						  pptr[i]->qmap,
						  pptr[i]->smap,
						  output_format,
						  out);
		++output_cnt;
		if (output_cnt == num_output) break;
        //if(pptr[0]->read_id==50){printf("output_cnt is%d\n",output_cnt);}
	}
}

int result_combine(int readcount, int filecount, char *workpath, char *outfile, char *fastaq, int main_argc, char* main_argv[])
{
	char path[1024], buffer[1024];
	sprintf(path, "%s/chrindex.txt", workpath);
	FILE* chr_idx_file = fopen(path, "r");
	if (!chr_idx_file) { fprintf(stderr, "failed to open file %s for reading.\n", path); abort(); }
	int num_chr = 0;
	while(fgets(buffer, 1024, chr_idx_file)) ++num_chr;
	--num_chr;
	fastaindexinfo* chr_idx = (fastaindexinfo*)malloc(sizeof(fastaindexinfo) * num_chr);
	fseek(chr_idx_file, 0L, SEEK_SET);
	int i, flag;
	for (i = 0; i < num_chr; ++i)
	{
		flag = fscanf(chr_idx_file, "%ld\t%s\t%ld\n", &chr_idx[i].chrstart, chr_idx[i].chrname, &chr_idx[i].chrsize);
		assert(flag == 3);
	}
	fclose(chr_idx_file);
	chr_idx_file = NULL;
	
	fprintf(stderr, "output file name: %s\n", outfile);
	FILE* out = fopen(outfile, "w");
	char* out_buffer = (char*)malloc(8192);
	setvbuf(out, out_buffer, _IOFBF, 8192);
	
	if (output_format == FMT_SAM)
	{
		print_sam_header(out);
		print_sam_references(chr_idx, num_chr, out);
		print_sam_program(main_argc, main_argv, out);
	}
    
	const int trsize = num_candidates + 6;
	TempResult* pptr[trsize];
	for (i = 0; i < trsize; ++i) pptr[i] = create_temp_result();
	int num_results = 0;
	TempResult* trslt = create_temp_result();
	char* trf_buffer = (char*)malloc(8192);
	for (i = 1; i <= filecount; ++i)
	{
		sprintf(path, "%s/%d.r", workpath, i);
		FILE* thread_results_file = fopen(path, "r");
		if (!thread_results_file) { fprintf(stderr, "failed to open file %s for reading.\n", path); abort(); }
		setvbuf(thread_results_file, trf_buffer, _IOFBF, 8192);
		num_results = 0;
		int rok = load_temp_result(trslt, thread_results_file);
		if (rok) { copy_temp_result(trslt, pptr[num_results]); ++num_results; }
		int last_qid = pptr[0]->read_id;
		while (rok)
		{
			rok = load_temp_result(trslt, thread_results_file);
			if (!rok) break;

			if (trslt->read_id != last_qid)
			{
				output_query_results(chr_idx, num_chr, pptr, num_results, out);
				num_results = 0;
			}
			
			last_qid = trslt->read_id;
			copy_temp_result(trslt, pptr[num_results]);
			++num_results;
		}
		if (num_results) output_query_results(chr_idx, num_chr, pptr, num_results, out);
		
		fclose(thread_results_file);
	}
	
	for (i = 0; i < trsize; ++i) pptr[i] = destroy_temp_result(pptr[i]);
	destroy_temp_result(trslt);
	fclose(out);
	free(out_buffer);
	free(trf_buffer);
	free(chr_idx);
	
	return 0;
}


extern int meap_ref_impl_large(int, int, int);


int delete_mini_result(TempResult **pptr,TempResult **out_pptr,int lenA, int lenB){
    int ndelcount=0;
    bool bsame=false;
    int j;
    for(int i=0;i<lenA;i++){
        bsame=false;
        for (int k=0;k<lenB;k++){
            if(pptr[i]==out_pptr[k]){
                ndelcount++;
                bsame=true;
                break;
            }
        }
        if(false==bsame){
            pptr[j++]==pptr[i];
        }
    }
    return lenA-ndelcount;
}
//**********这部分是将第一次改和第二次改的结果做一次整合****************
void polish_result(const char *workpath,int filecount,int refcount,char  *refoutfile,int main_argc, char* main_argv[]){
    char path[200],path2[200];FILE *thread_file; FILE **up_file;int num_count=0;char buffer[1024];
    char *trbuffer=(char *)malloc(8192);char tempstr[200];int temp1,temp2;int r_k;
    int num_results=0;char *ref_name;char *ref_name2;
    const int trsize=num_candidates + 6;
    int ref_trsize=refcount *5;
    TempResult *pptr[trsize];
    TempResult *out_pptr[trsize];
    //TempResult *Temp_pptr[trsize];
    for (int i = 0; i < trsize; ++i) {
        pptr[i] = create_temp_result();
        out_pptr[i] = create_temp_result();
    }
    TempResult *trslt=create_temp_result();
    char* trf_buffer = (char*)malloc(8192);
    
    
    
    TempResult *refpptr[big_size];
    FILE *thread_ref_file;int num_ref_results=0;
    TempResult *trslt1=create_temp_result();
    char* trf_buffer1 = (char*)malloc(8192);
    for (int i = 0; i < big_size; ++i) refpptr[i] = create_temp_result();
    for(int i=1;i<=filecount;++i){
        sprintf(path,"%s/ref%d.r",workpath,i);
        thread_ref_file=fopen(path,"r");
        setvbuf(thread_ref_file,trf_buffer1,_IOFBF,8192);
        int rok=load_temp_result(trslt1, thread_ref_file);
        if(rok){
            copy_temp_result(trslt1,refpptr[num_ref_results]);
            num_ref_results++;
        }
        
        while(rok){
            rok=load_temp_result(trslt1, thread_ref_file);
            if(!rok)break;
            copy_temp_result(trslt1,refpptr[num_ref_results]);
            ++num_ref_results;
        }
        fclose(thread_ref_file);
    }
    printf("num_ref_results is %d\n",num_ref_results);//到这里没问题
    int **result_database;int  pre_id,current_id;int kkk=0;int *point_arr;
    result_database=(int **)calloc(1000000,sizeof(int*));
    for(int uul=0;uul<1000000;uul++){
        result_database[uul]=(int*)calloc(40,sizeof(int));
    }
    pre_id=refpptr[0]->read_id;
    for(int i=0;i<num_ref_results;i++){
        current_id=refpptr[i]->read_id;
        if(pre_id==current_id){
            kkk++;
        }
        else{
            kkk=0;
            pre_id=current_id;
        }
        result_database[current_id][kkk]=i;
    }
    //printf("database is%d\n",result_database[10][2]);
    sprintf(path2, "%s/chrindex.txt", workpath);
    FILE* chr_idx_file = fopen(path2, "r");
    if (!chr_idx_file) { fprintf(stderr, "failed to open file %s for reading.\n", path); abort(); }
    int num_chr = 0;
    while(fgets(buffer, 1024, chr_idx_file)) ++num_chr;//参考基因组
    --num_chr;
    up_file=(FILE **)malloc(filecount*sizeof(FILE *));
    fastaindexinfo* chr_idx = (fastaindexinfo*)malloc(sizeof(fastaindexinfo) * num_chr);
    fseek(chr_idx_file, 0L, SEEK_SET);
    int i, flag,flag2=0;
    
    for (i = 0; i < num_chr; ++i)
    {
        flag = fscanf(chr_idx_file, "%ld\t%s\t%ld\n", &chr_idx[i].chrstart, chr_idx[i].chrname, &chr_idx[i].chrsize);
        assert(flag == 3);
    }
    fprintf(stderr, "output file name: %s\n", refoutfile);
    
    FILE* out = fopen(refoutfile, "w");
    if (output_format == FMT_SAM)
    {
        print_sam_header(out);
        print_sam_references(chr_idx, num_chr, out);
        print_sam_program(main_argc, main_argv, out);
    }
    char* out_buffer = (char*)malloc(8192);
    setvbuf(out, out_buffer, _IOFBF, 8192);
    fclose(chr_idx_file);
    chr_idx_file = NULL;
    for(int i=0;i<filecount;i++){
        sprintf(tempstr,"%s/up%d.r",workpath,filecount);
        //up_file[i]=fopen(tempstr,"w");
    }//初始化文件
    for(int ww=1;ww<=filecount;ww++)
    {
        sprintf(path,"%s/%d.r",workpath,ww);
        
        thread_file=fopen(path,"r");
        
        if(!thread_file){fprintf(stderr,"failed to open file %s for reading\n",path);abort();}
        
        setvbuf(thread_file,trbuffer,_IOFBF,8192);
        num_results=0;int judg=0;
        
        int rok=load_temp_result(trslt,thread_file);
        if(rok){copy_temp_result(trslt,pptr[num_results]);++num_results;}
        int last_id=pptr[0]->read_id;int formal_id;int formal_loc;int org_sta,org_end,org_ref_start,org_ref_end;int ref_sid2;int temp_sb,temp_se;
        int or_sta,or_end;int *vote;int *mark;int maxi=-1;int maxi_vote=-1;int p_num=0; int delete_flag=0;int ref_sid3;int flag3=0;int flag4=0;
        int ref_sid4;int sid;
        vote=(int *)malloc(16*sizeof(int));
        mark=(int *)malloc(16*sizeof(int));
        for(int i=0;i<16;i++){
            vote[i]=0;
            mark[i]=0;
        }
        //printf("test\n");
        int ref_size;
        while(rok){
            rok=load_temp_result(trslt, thread_file);
            if(!rok)break;
            if(trslt->read_id!=last_id){
                for(int i=0;i<num_results;i++){
                    for(int j=i+1;j<num_results;j++){
                        int sid=get_chr_id(chr_idx, num_chr, pptr[i]->sb);
                        ref_sid2=get_chr_id(chr_idx, num_chr, pptr[j]->sb);
                        if(sid==ref_sid2){
                            if (labs(pptr[j]->qb-pptr[i]->qb)>1000&&labs(pptr[i]->qe-pptr[j]->qe)>1000) {
                                if(pptr[i]->sb!=pptr[j]->sb){
                                    if (fabs((pptr[i]->qb-pptr[j]->qb)/(pptr[i]->sb-pptr[j]->sb)-1)<0.9) {
                                        vote[i]=vote[i]+1;
                                        vote[j]=vote[j]+1;
                                        mark[i]=1;
                                        mark[j]=1;
                                    }
                                }
                            }
                        }
                    }
                }
                
                for(int k=0;k<num_results;k++){
                    delete_flag=0;
                    
                    if(mark[k]==0){
                        out_pptr[p_num]=pptr[k];
                        ++p_num;
                    }
                    if (mark[k]==1){//&
                        maxi_vote=vote[k];
                        maxi=k;
                        sid=get_chr_id(chr_idx, num_chr, pptr[k]->sb);
                        for(int p=k+1;p<num_results;p++){//*
                            ref_sid2=get_chr_id(chr_idx, num_chr, pptr[p]->sb);
                            if(sid==ref_sid2){
                                if(labs(pptr[p]->qb-pptr[k]->qb)<1500){
                                    mark[p]=2;
                                    delete_flag=1;
                                    if(labs((pptr[k]->qe-pptr[k]->qb)-(pptr[p]->qe-pptr[p]->qb))>3000){
                                        if((pptr[k]->qe-pptr[k]->qb)>(pptr[p]->qe-pptr[p]->qb)){ maxi=k;}
                                        else{maxi=p;}
                                    }else{
                                        if(maxi_vote>vote[p]){
                                            maxi=k;
                                        }
                                        if (maxi_vote==vote[p]) {
                                            if((pptr[k]->qe-pptr[k]->qb)>(pptr[p]->qe-pptr[p]->qb)){
                                                maxi=k;
                                            }
                                            else{
                                                maxi=p;
                                            }
                                        }
                                        if(maxi_vote<vote[p]){
                                            maxi=p;
                                            maxi_vote=vote[p];
                                        }
                                    }
                                }
                                else{
                                    for(int q=p+1;q<num_results;q++){
                                        ref_sid3=get_chr_id(chr_idx, num_chr, pptr[q]->sb);
                                        
                                        if(sid == ref_sid3){
                                            flag3=1;//说明不是最后一个
                                            break;
                                        }
                                       
                                    }
                                    if(flag3==0){
                                        out_pptr[p_num]=pptr[k];
                                        ++p_num;
                                    }
                                    flag3=0;
                                }
                            }
                        }//*
                        if(delete_flag==1){
                            out_pptr[p_num]=pptr[maxi];
                            ++p_num;}
                        for(int w=k+1;w<num_results;w++){
                            
                            ref_sid4=get_chr_id(chr_idx, num_chr, pptr[w]->sb);
                            if(sid ==ref_sid4){
                                flag4=1;
                                break;//这也说明不是最后一个
                            }
                           
                        }
                        if(flag4==0){
                            out_pptr[p_num]=pptr[k];
                            ++p_num;
                        }
                        flag4=0;
                    }//&
                    
                }
                
                output_query_results(chr_idx, num_chr, out_pptr,p_num, out);//shuchu
                num_results=0;
                p_num=0;
                
               for(int i=0;i<16;i++){mark[i]=0;vote[i]=0;}
            }
            
            last_id = trslt->read_id;
            copy_temp_result(trslt, pptr[num_results]);
            ++num_results;
        }
            if (num_results) output_query_results(chr_idx, num_chr, pptr, num_results, out);
            fclose(thread_file);
    }
    for(int i=0;i<trsize;++i)pptr[i]=destroy_temp_result(pptr[i]);
    for(int i=0;i<big_size;++i)refpptr[i]=destroy_temp_result(refpptr[i]);
    destroy_temp_result(trslt);
    free(chr_idx);
    free(up_file);
    free(out_buffer);
    free(trf_buffer1);
    free(result_database);
    //free(refpptr);
    
}

#define __run_system(cmd) \
	do { \
	int __rc_status = system(cmd); \
	if (__rc_status != 0) { fprintf(stderr, "[%s, %u] system() error. Error code is %d.\n", __func__, __LINE__, __rc_status); exit(1); } \
} while (0);
int result_combine3(int readcount, int filecount, char *workpath, char *outfile, char *fastaq, int main_argc, char* main_argv[])
{
    char path[1024], buffer[1024];
    sprintf(path, "%s/chrindex.txt", workpath);
    FILE* chr_idx_file = fopen(path, "r");
    if (!chr_idx_file) { fprintf(stderr, "failed to open file %s for reading.\n", path); abort(); }
    int num_chr = 0;
    while(fgets(buffer, 1024, chr_idx_file)) ++num_chr;
    --num_chr;
    fastaindexinfo* chr_idx = (fastaindexinfo*)malloc(sizeof(fastaindexinfo) * num_chr);
    fseek(chr_idx_file, 0L, SEEK_SET);
    int i, flag;
    for (i = 0; i < num_chr; ++i)
    {
        flag = fscanf(chr_idx_file, "%ld\t%s\t%ld\n", &chr_idx[i].chrstart, chr_idx[i].chrname, &chr_idx[i].chrsize);
        assert(flag == 3);
    }
    fclose(chr_idx_file);
    chr_idx_file = NULL;
    
    fprintf(stderr, "output file name: %s\n", outfile);
    FILE* out = fopen(outfile, "w");
    char* out_buffer = (char*)malloc(8192);
    setvbuf(out, out_buffer, _IOFBF, 8192);
    
    if (output_format == FMT_SAM)
    {
        print_sam_header(out);
        print_sam_references(chr_idx, num_chr, out);
        print_sam_program(main_argc, main_argv, out);
    }
    
    const int trsize = num_candidates + 6;
    TempResult* pptr[trsize];
    for (i = 0; i < trsize; ++i) pptr[i] = create_temp_result();
    int num_results = 0;
    TempResult* trslt = create_temp_result();
    char* trf_buffer = (char*)malloc(8192);
    for (i = 1; i <= filecount; ++i)
    {
        sprintf(path, "%s/up%d.r", workpath, i);
        FILE* thread_results_file = fopen(path, "r");
        if (!thread_results_file) { fprintf(stderr, "failed to open reffile %s for reading.\n", path); abort(); }
        setvbuf(thread_results_file, trf_buffer, _IOFBF, 8192);
        num_results = 0;
        int rok = load_temp_result(trslt, thread_results_file);
        if (rok) { copy_temp_result(trslt, pptr[num_results]); ++num_results; }
        int last_qid = pptr[0]->read_id;
        while (rok)
        {
            rok = load_temp_result(trslt, thread_results_file);
            if (!rok) break;
            
            if (trslt->read_id != last_qid)
            {
                output_query_results(chr_idx, num_chr, pptr, num_results, out);
                num_results = 0;
            }
            
            last_qid = trslt->read_id;
            copy_temp_result(trslt, pptr[num_results]);
            ++num_results;
        }
        if (num_results) output_query_results(chr_idx, num_chr, pptr, num_results, out);
        
        fclose(thread_results_file);
    }
    
    for (i = 0; i < trsize; ++i) pptr[i] = destroy_temp_result(pptr[i]);
    destroy_temp_result(trslt);
    fclose(out);
    free(out_buffer);
    free(trf_buffer);
    free(chr_idx);
    return 0;
}

int main(int argc, char *argv[])
{
	prog_name = argv[0];
	
    char cmd[300], outfile[200],refoutfile[200];
    int corenum;
    long filelength;
    struct timeval tpstart, tpend;
    struct timeval mapstart, mapend;
    float timeuse;
    char saved[150], fastqfile[150], fastafile[150], tempstr[150],tempstr1[150];
    int  readcount;
    int  refcount;
    FILE *fid1, *fid2;

    gettimeofday(&tpstart, NULL);
    corenum = firsttask(argc, argv);
    printf("first task is sucess\n");//*******
    gettimeofday(&tpend, NULL);
    /*timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    fid1 = fopen("config.txt", "r");
	const char* read_results = NULL;
    read_results = fgets(saved, 150, fid1);
	assert(read_results);
    saved[strlen(saved) - 1] = '\0';
    read_results = fgets(fastafile, 150, fid1);
	assert(read_results);
    fastafile[strlen(fastafile) - 1] = '\0';
    read_results = fgets(fastqfile, 150, fid1);
	assert(read_results);
    fastqfile[strlen(fastqfile) - 1] = '\0';
    read_results = fgets(outfile, 150, fid1);
	assert(read_results);
    outfile[strlen(outfile) - 1] = '\0';
    read_results = fgets(refoutfile, 150, fid1);
    assert(read_results);
    refoutfile[strlen(refoutfile) - 1] = '\0';
    int num_read_items = fscanf(fid1, "%d %d %d\n", &corenum,&readcount,&refcount);
	assert(num_read_items == 3);
    fclose(fid1);
    filelength=get_file_size(fastafile);
    gettimeofday(&mapstart, NULL);
    printf("meap  pre is sucuess\n");
	meap_ref_impl_large(num_candidates, num_output, tech);
    gettimeofday(&mapend, NULL);
    timeuse = 1000000 * (mapend.tv_sec - mapstart.tv_sec) + mapend.tv_usec - mapstart.tv_usec;
    timeuse /= 1000000;
    sprintf(tempstr,"%s/0.fq",saved);
    sprintf(tempstr1,"%s/ref.fq",saved);
    result_combine(readcount, corenum, saved, outfile,tempstr, argc, argv);
    //result_combine2(refcount, corenum, saved, refoutfile,tempstr1, argc, argv);
    polish_result(saved,corenum,refcount,refoutfile,argc,argv);
    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    
    fid2 = fopen("config.txt", "a");
    fprintf(fid2, "The total Time : %f sec\n", timeuse);
    fclose(fid2);
    sprintf(cmd, "cp -r config.txt \"%s.config\"", outfile);
    sprintf(cmd, "cp -r config.txt \"%s.config\"", refoutfile);
    __run_system(cmd);*/
    
	
	return EXIT_SUCCESS;
}
