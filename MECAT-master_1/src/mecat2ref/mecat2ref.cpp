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
#define RM 100000

#include "output.h"
#include "../common/defs.h"

static const char* prog_name = NULL;
static const int kDefaultNumCandidates = 10;
static const int kDefaultNumOutput = 10;
static int num_candidates;
static int num_output;
static int output_format = FMT_REF;
static const int kDefaultOutputFormat = FMT_REF;
static int tech;
static const int kDefaultTech = TECH_PACBIO;
typedef struct{
    long start;//k-mer在参考基因组出现的开始位置
    int readID;
    int frequency;//出现次数
    char string[13];
    int next[13];//next数组。K_mer做比对的时候，
}FRE_Kmer;
typedef struct {
    int readid;
    char onedata[100];
}longreadinfo;

typedef struct
{
	const char* reads;
	const char* reference;
	const char* wrk_dir;
	const char* output;
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
	fprintf(stderr, "%s [-d reads] [-r reference] [-o output] [-w working dir] [-t threads]", prog_name);
	fprintf(stderr, "\n\n");
	fprintf(stderr, "options:\n");
	fprintf(stderr, "-d <string>\treads file name\n");
	fprintf(stderr, "-r <string>\treference file name\n");
	fprintf(stderr, "-o <string>\toutput file name\n");
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
	while((opt_char = getopt(argc, argv, "d:r:w:o:t:n:b:m:x:")) != -1)
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
//自己加的
int chang_refer_fastqfile(const char *fastaq, const char *fenfolder)
{
    FILE *fp, *ot;
    int kk = 0,read_len;
    char read_name[1000], onedata[1000000], buff1[1000], buff2[1000000], *fq, tempstr[200], *oq,ch;
    sprintf(tempstr, "%s/refer.fq", fenfolder);
    fp = fopen(fastaq, "r");
    ot = fopen(tempstr, "w");
    fq = (char *)malloc(1000000000);
    setvbuf(fp, fq, _IOFBF, 1000000000);
    oq = (char *)malloc(1000000000);
    setvbuf(ot, oq, _IOFBF, 1000000000);
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
int * getNextArray(char *ms)
{
    int *next = (int*)malloc(sizeof(int)*(strlen(ms)+1));
    if(strlen(ms) == 1)
    {
        *next = -1;
        return next;
    }
    next[0] = -1;
    next[1] = 0;
    int pos = 2;
    int cn  = 0;
    while(pos <= strlen(ms))
    {
        if(ms[pos-1] == ms[cn])
            next[pos++] = ++cn;  //1
        else if(cn>0)
            cn = next[cn];       //2
        else
            next[pos++] = 0;     //3
    }
    //for(int i=0; i!=strlen(ms)+1; i++)
    //  printf("next[%d] = %d\n", i, next[i]);
    return next;
}
int KMPmatch(char *ms, char *str)
{
    int *next = getNextArray(ms);
    int count = 0;
    int msl  = strlen(ms);
    int strl = strlen(str);
    if(msl>strl)
        return 0;
    int mspos = 0;
    int strpos = 0;
    while(strpos < strl)
    {
        if(ms[mspos] == str[strpos])
        {
            mspos++,strpos++;           //4
            if(mspos == msl)
                count++;
            if(strpos == strl)
                return count;
        }
        else{
            if(mspos == 0){
                strpos++;           //5
                if((strl-strpos) < msl)
                    return count;
            }
            else
                mspos = next[mspos];   //6
        }
    }
    return count;
}
static long get_file_size(char *path)
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
static void filename(char *name,char *path1,int read_number){//path2是long read 文件
    long length;
    length=get_file_size(name);
    printf("%ld\n",length);
    FILE* file_path=fopen(name,"r");
    if(file_path!=NULL){
        printf("sucess\n");
    }
    FRE_Kmer *Kmer_index;
    int BC=20;int K_mer_length=12;
    long K_mer_numbe=length/12;
    Kmer_index=(FRE_Kmer *)malloc(sizeof(FRE_Kmer)*K_mer_numbe);//给K_mer分配内存空间
    for(int QQ=0;QQ<=K_mer_numbe;QQ++){
        Kmer_index[QQ].frequency=0;
        Kmer_index[QQ].readID=0;
        Kmer_index[QQ].string[0]='\0';
        Kmer_index[QQ].start=0;
        
    }// 初始化
    
    if(Kmer_index!=NULL)printf(" K_MER success\n");
    long number=K_mer_numbe;
    
    for(int h=0;h<number;h++){
        char line[1000000];
        char str[10000000];
        int *l;//l表示read 的长度
        int ii=0;int jj=0;
        while(fgets(line,sizeof(line),file_path)){
            sscanf(line,"%*d %*d %s",str);
            //printf("%s\n",str);
            for(int KK=0;KK<=strlen(str);KK++){
                if(ii<12){
                    Kmer_index[jj].string[ii]=str[KK];
                    ii++;}
                else{
                    jj++;ii=0;
                    KK=KK-11;
                    Kmer_index[jj].string[ii]=str[KK];
                    ii++;
                }
                // printf("%c",Kmer_index[jj].string[ii]=str[KK]);
            }
            
        } //初始化FRE_kmer
        //算相似度,block的长度是1000，K_mer的长度是20个。reference//K_mer长度为12？？？
    }
    FILE *fp2;
    fp2=fopen(path1,"r");
    int count2=read_number;//这个可以从源文件中找；readcountnumber是已知的
    longreadinfo *info;
    char line1[rm];
    int h=0;
    info=(longreadinfo *)malloc(sizeof(longreadinfo)*count2);
    if(info!=NULL)printf("success  info\n");
    if(fp2!=NULL)printf("success22\n");
    while(fgets(line1,sizeof(line1),fp2)){
        
        char str1[RM];
        sscanf(line1,"%*d %*d %s",str1);
        
        for(int y=0;y<strlen(str1);y++){
            info[h].onedata[y]=str1[y];
        }
        
        info[h].readid=h;
        
        h++;
        
    }
    
    int count=0;int ff=0;//出现次数
    
    printf("test%s\n",info[1].onedata);
    printf("test%d\n",info[0].readid);
    for(int o=0;o<K_mer_numbe;o++){
        for(int oo=0;oo<count2;oo++){
            count=KMPmatch(Kmer_index[o].string,info[oo].onedata);
            printf("%s\n",Kmer_index[o].string);
            printf("%s\n",info[oo].onedata);
            
        }
        ff=ff+count;
        Kmer_index[o].frequency=ff;
        printf("count is%d\n",ff);
        
        
    }
    
    
}

int firsttask(int argc, char *argv[])
{
	meap_ref_options* options = (meap_ref_options*)malloc(sizeof(meap_ref_options));
	int flag = param_read_t(argc, argv, options);
	if (flag == -1) { print_usage(); exit(1); }
	
    int readcount = chang_fastqfile(options->reads, options->wrk_dir);
    int refer_readcount=chang_refer_fastqfile(options->reference,options->wrk_dir);
    char path_[100];char path2_[100];
    sprintf(path_, "%s/0.fq",options->wrk_dir);
    sprintf(path2_, "%s/refer.fq",options->wrk_dir);
    filename(path2_,path_,readcount);
	char kkkkk[1024];
    sprintf(kkkkk, "config.txt");
    FILE* fileout = fopen(kkkkk, "w");
    fprintf(fileout, "%s\n%s\n%s\n%s\n%d\t%d\n", options->wrk_dir, options->reference, options->reads, options->output, options->num_cores, readcount);
    fclose(fileout);
	int corenum = options->num_cores;
	num_candidates = options->num_candidates;
	num_output = options->num_output;
	output_format = options->output_format;
	tech = options->tech;
	free(options);
    return (corenum);
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

long get_file_size(const char *path)
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

extern int meap_ref_impl_large(int, int, int);

#define __run_system(cmd) \
	do { \
	int __rc_status = system(cmd); \
	if (__rc_status != 0) { fprintf(stderr, "[%s, %u] system() error. Error code is %d.\n", __func__, __LINE__, __rc_status); exit(1); } \
} while (0);

int main(int argc, char *argv[])
{
	prog_name = argv[0];
	
    char cmd[300], outfile[200];
    int corenum;
    long filelength;
    struct timeval tpstart, tpend;
    struct timeval mapstart, mapend;
    float timeuse;
    char saved[150], fastqfile[150], fastafile[150], tempstr[150];
    int  readcount;
    FILE *fid1, *fid2;

    gettimeofday(&tpstart, NULL);
    corenum = firsttask(argc, argv);
    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
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
    int num_read_items = fscanf(fid1, "%d %d\n", &corenum,&readcount);
	assert(num_read_items == 2);
    fclose(fid1);
    filelength=get_file_size(fastafile);
    gettimeofday(&mapstart, NULL);
	meap_ref_impl_large(num_candidates, num_output, tech);
    gettimeofday(&mapend, NULL);
    timeuse = 1000000 * (mapend.tv_sec - mapstart.tv_sec) + mapend.tv_usec - mapstart.tv_usec;
    timeuse /= 1000000;

    sprintf(tempstr,"%s/0.fq",saved);
    result_combine(readcount, corenum, saved, outfile,tempstr, argc, argv);
    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    fid2 = fopen("config.txt", "a");
    fprintf(fid2, "The total Time : %f sec\n", timeuse);
    fclose(fid2);
    sprintf(cmd, "cp -r config.txt \"%s.config\"", outfile);
    __run_system(cmd);
	
	return EXIT_SUCCESS;
}
