#include "mecat2ref_defs.h"
#include "output.h"
#include "mecat2ref_aux.h"
#include "../common/diff_gapalign.h"
#include "../common/xdrop_gapalign.h"

#include <algorithm>
using namespace std;

static int MAXC = 0;
static int TECH = TECH_PACBIO;
static int num_output = MAXC;
static const double ddfs_cutoff_pacbio = 0.25;
static const double ddfs_cutoff_nanopore = 0.1;
static double ddfs_cutoff = ddfs_cutoff_pacbio;

static pthread_t *thread;
static int threadnum=2;
static FILE **outfile;
static pthread_mutex_t mutilock; 
static int runnumber=0,runthreadnum=0, readcount,terminalnum;
static int *countin;
static long **databaseindex,*allloc,seqcount,sumcount;
static long **databaseindex1,*allloc1,seqcount1,sumcount1;
static char *REFSEQ;
static char *savework,workpath[300],fastqfile[300];
static int *countin1;
static int seed_len=13;
static int index_count=67108864;
static char *read_REFESQ;
static char *save_work;
static ReadFasta *readinfo;
int mavalue[2000];
int similarity_count;
int read_kmer;
float count_value;
typedef struct {
    char *read_string;
}read_info;
typedef struct{
    double simm;//相似度
    int k_count;//long_read出现的次数
    float LDF;//词频
    float TF;
    float vote;
    int r_count;//参考基因里出现的次数
} sim;
static sim *sc;
static read_info *info;

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
    return(sumval);
}

static int transnum_buchang(char *seqm,int *value,int *endn,int len_str,int readnum,int BC)
{
    int eit=0,temp;
    int i,j,start,num;
    num=(len_str-readnum)/BC+1;
    *endn=(len_str-readnum)%BC;
    //if((len_str%readnum)>0)num=num+1;
    for(i=0; i<num; i++)
    {
        eit=0;
        start=i*BC;//滑动窗口
        //if(i==num-1){eit=0;start=len_str-readnum;}
        for(j=0; j<readnum; j++)
        {
            temp=atcttrans(seqm[start+j]);
            if(temp==4)
            {
                eit=-1;
                break;
            }
            eit=eit<<2;
            eit=eit+temp;
        }
        value[i]=eit;
    }
    return(num);//extract_number
}

static void insert_loc(struct Back_List *spr,int loc,int seedn,float len)
{//insert_loc(temp_spr,u_k,k+1,BC)
    int list_loc[SI],list_score[SI],list_seed[SI],i,j,minval,mini;
    for(i=0; i<SM; i++)
    {
        list_loc[i]=spr->loczhi[i];
        list_seed[i]=spr->seedno[i];
        list_score[i]=0;
    }
    list_loc[SM]=loc;
    list_seed[SM]=seedn;//seednumber
    list_score[SM]=0;//SM 和SI 是20 和21
    mini=-1;
    minval=10000;
    for(i=0; i<SM; i++)for(j=i+1; j<SI; j++)if(list_seed[j]-list_seed[i]>0&&list_loc[j]-list_loc[i]>0&&fabs((list_loc[j]-list_loc[i])/((list_seed[j]-list_seed[i])*len)-1.0)<ddfs_cutoff)
            {
                list_score[i]++;
                list_score[j]++;
            }
    for(i=0; i<SI; i++)if(minval>list_score[i])
        {
            minval=list_score[i];
            mini=i;
        }//找出最低的分数
    if(minval==SM)
    {
        spr->loczhi[SM-1]=loc;
        spr->seedno[SM-1]=seedn;
    }
    else if(minval<SM&&mini<SM)
    {
        for(i=mini; i<SM; i++)
        {
            spr->loczhi[i]=list_loc[i+1];
            spr->seedno[i]=list_seed[i+1];
        }
        spr->score--;//删掉最低一个
    }
}
static void build_read_index(const char *path){
    unsigned int eit,temp;long start;
    char tempstr[200];
    sprintf(tempstr, "%s/0.fq",path);
    int leftnum;
    leftnum=34-2*seed_len;
    FILE  *fp;
    fp=fopen(tempstr,"r");
    char *seq;
    int length=get_file_size(path);
    printf(" read length is %d",length);
    
    read_REFESQ=(char *)malloc((1000000000+200000)*sizeof(char));
    seq=read_REFESQ;
    
    int templen;
    char *pre;int lenl=0;
    int flag;int readno,readlen;int read_count;
    pre=save_work;
    info=(read_info *)malloc((100000+2)*sizeof(read_info));
    int lenth_count=0;int read_len; int temp_len;
    int lenth2_count=0;
    while((flag=fscanf(fp,"%d\t%d\t%s\n",&readno,&readlen,pre))!=EOF&&read_count<100000&&lenl<1000000000)
    {
        
        info[read_count].read_string=pre;
        read_len=strlen(pre);
        lenl=lenl+read_len+1;
        pre=pre+read_len+1;
        read_count++;
        
        
    }
    for(int i=0;i<read_count;i++){
        temp_len=strlen(info[i].read_string);
        for(int j=0;j<temp_len;j++){
            seq[lenth2_count]=info[i].read_string[j];
            lenth2_count++;
        }
        
    }
    seq[lenth2_count+1]='\0';
    printf("%s",seq);
    int actual_len=strlen(seq);
    seqcount1=actual_len;
    seq[actual_len+1]='\0';
    //printf("Constructing look-up table...\n");
    countin1=(int *)malloc((index_count)*sizeof(int));
    for(int i=0; i<index_count; i++)countin1[i]=0;
    
    // Count the number
    eit=0;
    start=0;
    for(int i=0; i<seqcount1; i++)
    {
        //printf("%c",seq[i]);
        if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4)
        {
            eit=0;
            start=0;
            continue;
        }
        temp=atcttrans(seq[i]);
        if(start<seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            //printf("countin is %d\n",eit);
            
        }
        else if(start>=seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            countin1[eit]=countin1[eit]+1;
            //printf("countin is %d\n",eit);//存的是countin
            eit=eit<<leftnum;
            eit=eit>>leftnum;
            // printf("eit is %d\n",eit);
        }
        
    }
    
    
    //Max_index
    sumcount1=sumvalue_x(countin1,index_count);
    read_kmer=sumcount1;
    allloc1=(long *)malloc(sumcount1*sizeof(long));
    databaseindex1=(long **)malloc((index_count)*sizeof(long*));
    //allocate memory
    sumcount1=0;
    for( int i=0; i<index_count; i++)
    {
        if(countin1[i]>0)
        {
            databaseindex1[i]=allloc1+sumcount1;
            sumcount1=sumcount1+countin1[i];
            countin1[i]=0;
        }
        else databaseindex1[i]=NULL;
        
    }
    
    // printf("xiao");//10834098
    
    //constructing the look-up table
    eit=0;
    start=0;
    for(int i=0; i<seqcount1; i++)
    {
        //printf("%c\n",seq[i]);
        if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4)
        {
            eit=0;
            start=0;
            continue;
        }
        temp=atcttrans(seq[i]);
        if(start<seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            // printf("eit2%d\n",eit);
        }
        else if(start>=seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            
            if(databaseindex1[eit]!=NULL)
            {
                countin1[eit]=countin1[eit]+1;
                databaseindex1[eit][countin1[eit]-1]=i+2-seed_len;
                
            }
            eit=eit<<leftnum;
            eit=eit>>leftnum;
            printf("%d\n",eit);
            
        }
        
    }
    
    
    
    
}

static void creat_ref_index(char *fastafile)
{
    unsigned int eit,temp;
    int  indexcount=0,leftnum=0;
    long length,count,i,start, rsize = 0;
    FILE *fasta,*fastaindex;
    char *seq,ch,nameall[200];
    if(seed_len==14)indexcount=268435456;
    else if(seed_len==13)indexcount=67108864;//2的26次方
    else if(seed_len==12)indexcount=16777216;
    else if(seed_len==11)indexcount=4194304;
    else if(seed_len==10)indexcount=1048576;
    else if(seed_len==9)indexcount=262144;
    else if(seed_len==8)indexcount=65536;
    else if(seed_len==7)indexcount=16384;
    else if(seed_len==6)indexcount=4096;
    leftnum=34-2*seed_len;
    //read reference seq
    length=get_file_size(fastafile);
    fasta=fopen(fastafile, "r");
    sprintf(nameall,"%s/chrindex.txt",workpath);
    fastaindex=fopen(nameall,"w");

    REFSEQ=(char *)malloc((length+1000)*sizeof(char));
    seq=REFSEQ;
    for (ch=getc(fasta),count=0; ch!=EOF; ch=getc(fasta))
    {
        if(ch=='>')
        {
            assert(fscanf(fasta,"%[^\n]s",nameall) == 1);
			if (rsize) fprintf(fastaindex, "%ld\n", rsize);
			rsize = 0;
			for(i=0;i<strlen(nameall);i++)if(nameall[i]==' '||nameall[i]=='\t')break;
			nameall[i]='\0';
            fprintf(fastaindex,"%ld\t%s\t",count,nameall);
        }
        else if(ch!='\n'&&ch!='\r')
        {
            if(ch>'Z')ch=toupper(ch);
            seq[count]=ch;
            count=count+1;
			++rsize;
        }
    }
    fclose(fasta);
	fprintf(fastaindex, "%ld\n", rsize);
    fprintf(fastaindex,"%ld\t%s\n",count,"FileEnd");
    seq[count]='\0';
    fclose(fastaindex);
    seqcount=count;
    similarity_count=(seqcount-12)/200+1;
    sc=(sim *)malloc(similarity_count*sizeof(sim));
    for(int k=0;k<similarity_count;k++){
        sc[k].k_count=0;
        sc[k].simm=0;
        sc[k].TF=0;
        sc[k].LDF=0;//词频
        sc[k].r_count=0;
        sc[k].vote=0;
    }//初始化
    printf("%ld\n",seqcount);
//printf("Constructing look-up table...\n");
    countin=(int *)malloc((indexcount)*sizeof(int));
    for(i=0; i<indexcount; i++)countin[i]=0;

// Count the number
    eit=0;
    start=0;
    for(i=0; i<seqcount; i++)
    {
        if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4)
        {
            eit=0;
            start=0;
            continue;
        }
        temp=atcttrans(seq[i]);
        if(start<seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            //printf("%d\n",eit);
        }
        else if(start>=seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            countin[eit]=countin[eit]+1;
            eit=eit<<leftnum;
            eit=eit>>leftnum;
            // printf("this is eit%d\n",eit);
        }
       
        //printf("%c",seq[i]);
    }

  int nn;//表示区域
//Max_index
    sumcount=sumvalue_x(countin,indexcount);//有效Index K_mer的数量
    count_value=sumcount;
    allloc=(long *)malloc(sumcount*sizeof(long));
    databaseindex=(long **)malloc((indexcount)*sizeof(long*));
//allocate memory
    sumcount=0;
    for(i=0; i<indexcount; i++)
    {
        if(countin[i]>0)
        {
            databaseindex[i]=allloc+sumcount;
            sumcount=sumcount+countin[i];
            countin[i]=0;
        }
        else databaseindex[i]=NULL;
       
    }

    // printf("xiao");//10834098

//constructing the look-up table
    eit=0;
    start=0;
    for(i=0; i<seqcount; i++)
    {
        //printf("%c",seq[i]);
        if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4)
        {
            eit=0;
            start=0;
            continue;
        }
        temp=atcttrans(seq[i]);
        if(start<seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            //printf("eit2%d\n",eit);
        }
        else if(start>=seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;

            if(databaseindex[eit]!=NULL)
            {
                countin[eit]=countin[eit]+1;
                databaseindex[eit][countin[eit]-1]=i+2-seed_len;//存的位置
            }
            nn=(i-12)/200+1;//按照200划分，
            if(countin1[eit]>0){
                sc[nn].k_count=sc[nn].k_count+countin1[eit];//在long_read里面出现的次数
            };
            
            eit=eit<<leftnum;
            eit=eit>>leftnum;
            //printf(" thiS IS eit2%d\n",eit);
        }
        
    }
}
static void get_vote(){
    int eit=0;
    int temp=0;
    int i=0;char *seq;sim *sc1;char *readseq;
    readseq=read_REFESQ;
    seq=REFSEQ;
    sc1=sc;
    int start=0;//num 有关
    int leftnum=8;int nn=0;
    
    
    for(int j=0;j<similarity_count;j++){
        
        if(sc1[j].k_count>0){
            sc1[j].LDF=log((read_kmer)/sc1[j].k_count);
            printf("LDF is %f\n",sc1[j].LDF);}
    }//在Longread里面出现的次数
    
    for(i=0; i<seqcount; i++)
    {
        
        if(seq[i]=='N'||(temp=atcttrans(seq[i]))==4)
        {
            eit=0;
            start=0;
            continue;
        }
        temp=atcttrans(seq[i]);
        if(start<seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
        }
        else if(start>=seed_len-1)
        {
            eit=eit<<2;
            eit=eit+temp;
            start=start+1;
            nn=(i-12)/200+1;
            if(countin[eit]>0){
                printf("count in is %d\n",countin[eit]);
                sc1[nn].r_count=sc1[nn].r_count+countin[eit];
                
            }//在参考基因里出现的次数
            
            eit=eit<<leftnum;
            eit=eit>>leftnum;
        }
        
    }
    for(int j=0;j<similarity_count;j++){
        if(count_value>0){
            printf("r_count is %d\n",sc1[j].r_count);
            sc1[j].TF=sc1[j].r_count/count_value;
            printf("TF is %f\n",sc1[j].TF);}
        sc1[j].vote=sc1[j].TF*sc1[j].LDF;
        printf("vote is %f\n",sc1[j].vote);
        
    }
}



static void reference_mapping(int threadint)
{
    int cleave_num,read_len;
    int mvalue[20000],flag_end;
    long *leadarray,u_k,s_k,loc;
    int count1=0,i,j,k,templong,read_name;
    struct Back_List *database,*temp_spr,*temp_spr1;
    int repeat_loc = 0,*index_list,*index_spr;
    long location_loc[4],left_length1,right_length1,left_length2,right_length2,loc_list,start_loc;
    short int *index_score,*index_ss;
    int temp_list[200],temp_seedn[200],temp_score[200];
    int localnum,read_i,read_end,fileid;
    int endnum,ii;
    char *seq,*onedata,onedata1[RM],onedata2[RM],FR;
    int cc1,canidatenum,loc_seed;
    int num1,num2,BC;
    int low,high,mid,seedcount;
    candidate_save canidate_loc[MAXC],canidate_temp;
    seq=REFSEQ;
    j=seqcount/ZV+5;
	
	int* fwd_index_list = (int*)malloc(sizeof(int) * j);
	short* fwd_index_score = (short*)malloc(sizeof(short) * j);
	Back_List* fwd_database = (Back_List*)malloc(sizeof(Back_List) * j);
	int* rev_index_list = (int*)malloc(sizeof(int) * j);
	short* rev_index_score = (short*)malloc(sizeof(short) * j);
	Back_List* rev_database = (Back_List*)malloc(sizeof(Back_List) * j);//和参考基因有关的
	for (i = 0; i < j; ++i) {
		fwd_database[i].score = 0;
		fwd_database[i].score2 = 0;
		fwd_database[i].index = -1;
		rev_database[i].score = 0;
		rev_database[i].score2 = 0;
		rev_database[i].index = -1;
	}
	int fnblk, rnblk;
	int* pnblk;
	AlignInfo alns[MAXC + 6];
	int naln;
	TempResult results[MAXC + 6];
	int nresults;
	long aln_bytes = 1;
	aln_bytes = aln_bytes * 2 * (MAXC + 6) * MAX_SEQ_SIZE;
	char* aln_seqs = new char[aln_bytes];
	u_k = 0;
	for (int i = 0; i < MAXC + 6; ++i) {
		results[i].qmap = aln_seqs + u_k;
		u_k += MAX_SEQ_SIZE;
		results[i].smap = aln_seqs + u_k;
		u_k += MAX_SEQ_SIZE;
	}
	
	vector<char> qstr;
	vector<char> tstr;
	GapAligner* aligner = NULL;
	if (TECH == TECH_PACBIO) {
		aligner = new DiffAligner(0);
	} else if (TECH == TECH_NANOPORE) {
		aligner = new XdropAligner(0);
	} else {
		ERROR("TECH must be either %d or %d", TECH_PACBIO, TECH_NANOPORE);
	}

    fileid=1;
    while(fileid)
    {
        pthread_mutex_lock(&mutilock);
        localnum=runnumber;
        runnumber++;
        pthread_mutex_unlock(&mutilock);
        if(localnum>=terminalnum)
        {
            fileid=0;
            break;
        }
        if(localnum==terminalnum-1)read_end=readcount;
        else read_end=(localnum+1)*PLL;
        for(read_i=localnum*PLL; read_i<read_end; read_i++)
        {
            read_name=readinfo[read_i].readno;
            read_len=readinfo[read_i].readlen;
            strcpy(onedata1,readinfo[read_i].seqloc);

            canidatenum=0;
            for(ii=1; ii<=2; ii++)
            {
                read_len=strlen(onedata1);
                BC=5+(read_len/1000);
                if(BC>20)BC=20;
                if(ii==1) {
					onedata=onedata1;
					index_list = fwd_index_list;// int repeat_loc = 0,*index_list,*index_spr;
					index_score = fwd_index_score;//short int *index_score,*index_ss;
					database = fwd_database;
					pnblk = &fnblk;
				} else if(ii==2){
					index_list = rev_index_list;
					index_score = rev_index_score;
					database = rev_database;
					pnblk = &rnblk;
                    strcpy(onedata2,onedata1);
                    onedata=onedata2;
                    for(j=read_len-1,i=0; j>i; j--,i++)
                    {
                        FR=onedata[i];
                        onedata[i]=onedata[j];
                        onedata[j]=FR;
                    }
                    for(i=0; i<read_len; i++)
                    {
                        FR=onedata[i];
                        switch(FR)
                        {
                        case 'A':
                        {
                            onedata[i]='T';
                            break;
                        }
                        case 'T':
                        {
                            onedata[i]='A';
                            break;
                        }
                        case 'C':
                        {
                            onedata[i]='G';
                            break;
                        }
                        case 'G':
                        {
                            onedata[i]='C';
                            break;
                        }
                        }
                    }//将字符串掉个
                }
                endnum=0;
                read_len=strlen(onedata);
                cleave_num=transnum_buchang(onedata,mvalue,&endnum,read_len,seed_len,BC);
                j=0;//j=0
                index_spr=index_list;
                index_ss=index_score;
                endnum=0;
                for(k=0; k<cleave_num; k++)if(mvalue[k]>=0)
                    {
                        count1=countin[mvalue[k]];
                        //if(count1>20)continue;
                        leadarray=databaseindex[mvalue[k]];//leadarry存的位置
                        for(i=0; i<count1; i++,leadarray++)
                        {
                            templong=(*leadarray)/ZV;//位置信息，第几个block
                            u_k=(*leadarray)%ZV;//zv为1000，U_K是偏移量
                            if(templong>=0)
                            {
                                temp_spr=database+templong;//temp_spr的类型是back_list。   struct Back_List *database
                                if(temp_spr->score==0||temp_spr->seednum<k+1)
                                {
                                    loc=++(temp_spr->score);
                                    if(loc<=SM)
                                    {
                                        temp_spr->loczhi[loc-1]=u_k;//位置，block位置。
                                        temp_spr->seedno[loc-1]=k+1;
                                    }
                                    else insert_loc(temp_spr,u_k,k+1,BC);//删除分数最小的。保持在20个左右//(struct Back_List *spr,int loc,int seedn,float len)
                                    if(templong>0)s_k=temp_spr->score+(temp_spr-1)->score;
                                    else s_k=temp_spr->score;
                                    if(endnum<s_k)endnum=s_k;
                                    if(temp_spr->index==-1)
                                    {
                                        *(index_spr++)=templong;//int *index_list,*index_spr;
                                        *(index_ss++)=s_k;
                                        temp_spr->index=j;
                                        j++;//
                                    } else index_score[temp_spr->index]=s_k;
									temp_spr->score2 = temp_spr->score;
                                }
                                temp_spr->seednum=k+1;
                            }
                        }
                    }
				*pnblk = j;
                cc1=j;//index 数
                for(i=0,index_spr=index_list,index_ss=index_score; i<cc1; i++,index_spr++,index_ss++)if(*index_ss>6)//short int *index_score,*index_ss;
                    {
                        temp_spr=database+*index_spr;
                        if(temp_spr->score==0)continue;
                        s_k=temp_spr->score;
                        if(*index_spr>0)loc=(temp_spr-1)->score;
                        else loc=0;
                        start_loc=(*index_spr)*ZVL;
                        if(*index_spr>0)
                        {
                            loc=(temp_spr-1)->score;
                            if(loc>0)start_loc=(*index_spr-1)*ZVL;//记录位置的。
                        }
                        else loc=0;
                        if(loc==0)for(j=0,u_k=0; j<s_k&&j<SM; j++)
                            {
                                temp_list[u_k]=temp_spr->loczhi[j];
                                temp_seedn[u_k]=temp_spr->seedno[j];
                                u_k++;
                            }
                        else
                        {
                            k=loc;
                            u_k=0;
                            temp_spr1=temp_spr-1;
                            for(j=0; j<k&&j<SM; j++)
                            {
                                temp_list[u_k]=temp_spr1->loczhi[j];//int temp_list,temp_seedn类型
                                temp_seedn[u_k]=temp_spr1->seedno[j];
                                u_k++;
                            }
                            for(j=0; j<s_k&&j<SM; j++)
                            {
                                temp_list[u_k]=temp_spr->loczhi[j]+ZV;
                                temp_seedn[u_k]=temp_spr->seedno[j];
                                u_k++;
                            }
                        }
                        flag_end=find_location(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff);
                        if(flag_end==0)continue;
                        if(temp_score[repeat_loc]<6)continue;
                        canidate_temp.score=temp_score[repeat_loc];
                        loc_seed=temp_seedn[repeat_loc];
                        //loc_list=temp_list[repeat_loc];
                        location_loc[0]=start_loc+location_loc[0];
                        location_loc[1]=(location_loc[1]-1)*BC;
                        loc_list=location_loc[0];
                        left_length1=location_loc[0]+seed_len-1;
                        right_length1=seqcount-location_loc[0];
                        left_length2=location_loc[1]+seed_len-1;
                        right_length2=read_len-location_loc[1];
                        if(left_length1>=left_length2)num1=left_length2;
                        else num1=left_length1;
                        if(right_length1>=right_length2)num2=right_length2;
                        else num2=right_length1;
                        seedcount=0;
                        canidate_temp.loc1=location_loc[0];
                        canidate_temp.num1=num1;
                        canidate_temp.loc2=location_loc[1];
                        canidate_temp.num2=num2;
                        canidate_temp.left1=left_length1;
                        canidate_temp.left2=left_length2;
                        canidate_temp.right1=right_length1;
                        canidate_temp.right2=right_length2;
                        //find all left seed
                        for(u_k=*index_spr-2,k=num1/ZV,temp_spr1=temp_spr-2; u_k>=0&&k>=0; temp_spr1--,k--,u_k--)if(temp_spr1->score>0)
                            {
                                start_loc=u_k*ZVL;
								int scnt = min((int)temp_spr1->score, SM);
                                for(j=0,s_k=0; j < scnt; j++)if(fabs((loc_list-start_loc-temp_spr1->loczhi[j])/((loc_seed-temp_spr1->seedno[j])*BC*1.0)-1.0)<ddfs_cutoff)
                                    {
                                        seedcount++;
                                        s_k++;
                                    }
                                if(s_k*1.0/scnt>0.4)temp_spr1->score=0;
                            }
                        //find all right seed
                        for(u_k=*index_spr+1,k=num2/ZV,temp_spr1=temp_spr+1; k>0; temp_spr1++,k--,u_k++)if(temp_spr1->score>0)
                            {
                                start_loc=u_k*ZVL;
								int scnt = min((int)temp_spr1->score, SM);
                                for(j=0,s_k=0; j < scnt; j++)if(fabs((start_loc+temp_spr1->loczhi[j]-loc_list)/((temp_spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<ddfs_cutoff)
                                    {
                                        seedcount++;
                                        s_k++;
                                    }
                                if(s_k*1.0/scnt>0.4)temp_spr1->score=0;
                            }
                        canidate_temp.score=canidate_temp.score+seedcount;
                        if(ii==1)canidate_temp.chain='F';
                        else canidate_temp.chain='R';
                        //insert canidate position or delete this position
                        low=0;
                        high=canidatenum-1;
                        while(low<=high)
                        {
                            mid=(low+high)/2;
                            if(mid>=canidatenum||canidate_loc[mid].score<canidate_temp.score)high=mid-1;
                            else low=mid+1;
                        }
                        if(canidatenum<MAXC)for(u_k=canidatenum-1; u_k>high; u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
                        else for(u_k=canidatenum-2; u_k>high; u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
                        if(high+1<MAXC)canidate_loc[high+1]=canidate_temp;
                        if(canidatenum<MAXC)canidatenum++;
                        else canidatenum=MAXC;
                    }
            }

			naln = 0;
			nresults = 0;
            for(i=0; i<canidatenum; i++)
            {
				extend_candidate(canidate_loc[i], 
								 aligner, 
								 seq, 
								 seqcount,
								 onedata1, 
								 onedata2, 
								 qstr, 
								 tstr, 
								 read_name, 
								 read_len, 
								 alns, 
								 &naln,
								 results,
								 nresults);
            }
			
			rescue_clipped_align(alns, 
								  naln, 
								  results,
								  nresults,
								  aligner, 
								  seq, 
								  seqcount,
								  onedata1, 
								  onedata2, 
								  qstr, 
								  tstr, 
								  read_name, 
								  read_len, 
								  ZV, 
								  BC, 
								  fwd_database, 
								  rev_database,
								  ddfs_cutoff);
			
			output_results(alns, naln, results, nresults, num_output, outfile[threadint]);
			
			for (int t = 0; t < fnblk; ++t) {
				int bid = fwd_index_list[t];
				fwd_database[bid].score = 0;
				fwd_database[bid].score2 = 0;
				fwd_database[bid].index = -1;
			}
			for (int t = 0; t < rnblk; ++t) {
				int bid = rev_index_list[t];
				rev_database[bid].score = 0;
				rev_database[bid].score2 = 0;
				rev_database[bid].index = -1;
			}

			if (naln == 0)
            {
                canidatenum=0;
                for(ii=1; ii<=2; ii++)
                {
                    read_len=strlen(onedata1);
                    BC=5;
                    if(ii==1) {
						onedata=onedata1;
						index_list = fwd_index_list;
						index_score = fwd_index_score;
						database = fwd_database;
						pnblk = &fnblk;
					} else if(ii==2) {
						index_list = rev_index_list;
						index_score = rev_index_score;
						database = rev_database;
						pnblk = &rnblk;
                        strcpy(onedata2,onedata1);
                        onedata=onedata2;
                        for(j=read_len-1,i=0; j>i; j--,i++)
                        {
                            FR=onedata[i];
                            onedata[i]=onedata[j];
                            onedata[j]=FR;
                        }
                        for(i=0; i<read_len; i++)
                        {
                            FR=onedata[i];
                            switch(FR)
                            {
                            case 'A':
                            {
                                onedata[i]='T';
                                break;
                            }
                            case 'T':
                            {
                                onedata[i]='A';
                                break;
                            }
                            case 'C':
                            {
                                onedata[i]='G';
                                break;
                            }
                            case 'G':
                            {
                                onedata[i]='C';
                                break;
                            }
                            }
                        }
                    }

                    endnum=0;
                    read_len=strlen(onedata);
                    cleave_num=transnum_buchang(onedata,mvalue,&endnum,read_len,seed_len,BC);
                    j=0;
                    index_spr=index_list;
                    index_ss=index_score;
                    endnum=0;
                    for(k=0; k<cleave_num; k++)if(mvalue[k]>=0)
                        {
                            count1=countin[mvalue[k]];
                            //if(count1>20)continue;
                            leadarray=databaseindex[mvalue[k]];
                            for(i=0; i<count1; i++,leadarray++)
                            {
                                templong=(*leadarray)/ZVS;
                                u_k=(*leadarray)%ZVS;
                                if(templong>=0)
                                {
                                    temp_spr=database+templong;
                                    if(temp_spr->score==0||temp_spr->seednum<k+1)
                                    {
                                        loc=++(temp_spr->score);
                                        if(loc<=SM)
                                        {
                                            temp_spr->loczhi[loc-1]=u_k;
                                            temp_spr->seedno[loc-1]=k+1;
                                        }
                                        else insert_loc(temp_spr,u_k,k+1,BC);
                                        if(templong>0)s_k=temp_spr->score+(temp_spr-1)->score;
                                        else s_k=temp_spr->score;
                                        if(endnum<s_k)endnum=s_k;
                                        if(temp_spr->index==-1)
                                        {
                                            *(index_spr++)=templong;
                                            *(index_ss++)=s_k;
                                            temp_spr->index=j;
                                            j++;
                                        } else index_score[temp_spr->index]=s_k;
										temp_spr->score2 = temp_spr->score;
                                    }
                                    temp_spr->seednum=k+1;
                                }
                            }
                        }
					*pnblk = j;
                    cc1=j;
                    for(i=0,index_spr=index_list,index_ss=index_score; i<cc1; i++,index_spr++,index_ss++)if(*index_ss>4)
                        {
                            temp_spr=database+*index_spr;
                            if(temp_spr->score==0)continue;
                            s_k=temp_spr->score;
                            if(*index_spr>0)loc=(temp_spr-1)->score;
                            else loc=0;
                            start_loc=(*index_spr)*ZVSL;
                            if(*index_spr>0)
                            {
                                loc=(temp_spr-1)->score;
                                if(loc>0)start_loc=(*index_spr-1)*ZVSL;
                            }
                            else loc=0;
                            if(loc==0)for(j=0,u_k=0; j<s_k&&j<SM; j++)
                                {
                                    temp_list[u_k]=temp_spr->loczhi[j];
                                    temp_seedn[u_k]=temp_spr->seedno[j];
                                    u_k++;
                                }
                            else
                            {
                                k=loc;
                                u_k=0;
                                temp_spr1=temp_spr-1;
                                for(j=0; j<k&&j<SM; j++)
                                {
                                    temp_list[u_k]=temp_spr1->loczhi[j];
                                    temp_seedn[u_k]=temp_spr1->seedno[j];
                                    u_k++;
                                }
                                for(j=0; j<s_k&&j<SM; j++)
                                {
                                    temp_list[u_k]=temp_spr->loczhi[j]+ZVS;
                                    temp_seedn[u_k]=temp_spr->seedno[j];
                                    u_k++;
                                }
                            }
                            flag_end=find_location(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff);
                            if(flag_end==0)continue;
                            if(temp_score[repeat_loc]<6)continue;
                            canidate_temp.score=temp_score[repeat_loc];
                            loc_seed=temp_seedn[repeat_loc];
                            //loc_list=temp_list[repeat_loc];
                            location_loc[0]=start_loc+location_loc[0];
                            location_loc[1]=(location_loc[1]-1)*BC;
                            loc_list=location_loc[0];
                            left_length1=location_loc[0]+seed_len-1;
                            right_length1=seqcount-location_loc[0];
                            left_length2=location_loc[1]+seed_len-1;
                            right_length2=read_len-location_loc[1];
                            if(left_length1>=left_length2)num1=left_length2;
                            else num1=left_length1;
                            if(right_length1>=right_length2)num2=right_length2;
                            else num2=right_length1;
                            seedcount=0;
                            canidate_temp.loc1=location_loc[0];
                            canidate_temp.num1=num1;
                            canidate_temp.loc2=location_loc[1];
                            canidate_temp.num2=num2;
                            canidate_temp.left1=left_length1;
                            canidate_temp.left2=left_length2;
                            canidate_temp.right1=right_length1;
                            canidate_temp.right2=right_length2;
                            //find all left seed
                            for(u_k=*index_spr-2,k=num1/ZVS,temp_spr1=temp_spr-2; u_k>=0&&k>=0; temp_spr1--,k--,u_k--)if(temp_spr1->score>0)
                                {
                                    start_loc=u_k*ZVSL;
									int scnt = min((int)temp_spr1->score, SM);
                                    for(j=0,s_k=0; j < scnt; j++)if(fabs((loc_list-start_loc-temp_spr1->loczhi[j])/((loc_seed-temp_spr1->seedno[j])*BC*1.0)-1.0)<ddfs_cutoff)
                                        {
                                            seedcount++;
                                            s_k++;
                                        }
                                    if(s_k*1.0/scnt>0.4)temp_spr1->score=0;
                                }
                            //find all right seed
                            for(u_k=*index_spr+1,k=num2/ZVS,temp_spr1=temp_spr+1; k>0; temp_spr1++,k--,u_k++)if(temp_spr1->score>0)
                                {
                                    start_loc=u_k*ZVSL;
									int scnt = min((int)temp_spr1->score, SM);
                                    for(j=0,s_k=0; j < scnt; j++)if(fabs((start_loc+temp_spr1->loczhi[j]-loc_list)/((temp_spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<ddfs_cutoff)
                                        {
                                            seedcount++;
                                            s_k++;
                                        }
                                    if(s_k*1.0/scnt>0.4)temp_spr1->score=0;
                                }
                            canidate_temp.score=canidate_temp.score+seedcount;
                            if(ii==1)canidate_temp.chain='F';
                            else canidate_temp.chain='R';
                            //insert canidate position or delete this position
                            low=0;
                            high=canidatenum-1;
                            while(low<=high)
                            {
                                mid=(low+high)/2;
                                if(mid>=canidatenum||canidate_loc[mid].score<canidate_temp.score)high=mid-1;
                                else low=mid+1;
                            }
                            if(canidatenum<MAXC)for(u_k=canidatenum-1; u_k>high; u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
                            else for(u_k=canidatenum-2; u_k>high; u_k--)canidate_loc[u_k+1]=canidate_loc[u_k];
                            if(high+1<MAXC)canidate_loc[high+1]=canidate_temp;
                            if(canidatenum<MAXC)canidatenum++;
                            else canidatenum=MAXC;
                        }
                }

				naln = 0;
				nresults = 0;
                for(i=0; i<canidatenum; i++)
                {
					extend_candidate(canidate_loc[i], 
									 aligner, 
									 seq, 
									 seqcount,
									 onedata1, 
									 onedata2, 
									 qstr, 
									 tstr, 
									 read_name, 
									 read_len, 
									 alns, 
									 &naln, 
									 results,
									 nresults);
                }
				
				rescue_clipped_align(alns, 
									  naln, 
									  results,
									  nresults,
									  aligner, 
									  seq, 
									  seqcount,
									  onedata1, 
									  onedata2, 
									  qstr, 
									  tstr, 
									  read_name, 
									  read_len, 
									  ZVS, 
									  BC, 
									  fwd_database, 
									  rev_database,
									  ddfs_cutoff);
				
				output_results(alns, naln, results, nresults, num_output, outfile[threadint]);
				
				for (int t = 0; t < fnblk; ++t) {
					int bid = fwd_index_list[t];
					fwd_database[bid].score = 0;
					fwd_database[bid].score2 = 0;
					fwd_database[bid].index = -1;
				}
				for (int t = 0; t < rnblk; ++t) {
					int bid = rev_index_list[t];
					rev_database[bid].score = 0;
					rev_database[bid].score2 = 0;
					rev_database[bid].index = -1;
				}
            }
        }
    }
	delete aligner;
    free(fwd_database);
    free(fwd_index_list);
    free(fwd_index_score);
	free(rev_database);
    free(rev_index_list);
    free(rev_index_score);
	delete[] aln_seqs;
}


static void* multithread(void* arg)
{
    int localthreadno;
    pthread_mutex_lock(&mutilock);
    localthreadno=runthreadnum;
    runthreadnum++;
    pthread_mutex_unlock(&mutilock);
    reference_mapping(localthreadno);
	return NULL;
}

static int load_fastq(FILE *fq)
{
    int readlen,readno,sum=0,flag;
    char *pre;
    readcount=0;
    pre=savework;
    while((flag=fscanf(fq,"%d\t%d\t%s\n",&readno,&readlen,pre))!=EOF&&readcount<SVM&&sum<MAXSTR)
    {
        readinfo[readcount].seqloc=pre;
        readinfo[readcount].readno=readno;
        readlen=strlen(pre);
        readinfo[readcount].readlen=readlen;
        sum=sum+readlen+1;
        pre=pre+readlen+1;
        readcount++;
    }
    if(flag!=EOF)
        
    {
        readinfo[readcount].seqloc=pre;
        readinfo[readcount].readno=readno;
        readlen=strlen(pre);
        readinfo[readcount].readlen=readlen;
        readcount++;
        return(1);
    }
    else return(0);
}


int meap_ref_impl_large(int maxc, int noutput, int tech)
{
	MAXC = maxc;
	TECH = tech;
	num_output = noutput;
    char tempstr[300],fastafile[300];
    int corenum,readall;//readall是readcount 的条数
    int fileflag,threadno,threadflag;
    FILE *fp,*fastq;
    struct timeval tpstart, tpend;
    float timeuse;
    save_work=(char *)malloc((MAXSTR+RM)*sizeof(char));
    fp=fopen("config.txt","r");
    assert(fscanf(fp,"%s\n%s\n%s\n%s\n%d %d\n",workpath,fastafile,fastqfile,tempstr,&corenum,&readall) == 6);
    fclose(fp);
    threadnum=corenum;
    //building reference index
    gettimeofday(&tpstart, NULL);
    seed_len=13;
    build_read_index(workpath);
    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    fp = fopen("config.txt", "a");
    fprintf(fp, "The Building read Index Time: %f sec\n", timeuse);
    fclose(fp);
    //build read_long index
    gettimeofday(&tpstart, NULL);
    creat_ref_index(fastafile);
    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    fp = fopen("config.txt", "a");
    fprintf(fp, "The Building  Reference  Index Time: %f sec\n", timeuse);
    fclose(fp);
    get_vote();
    gettimeofday(&tpstart, NULL);

    savework=(char *)malloc((MAXSTR+RM)*sizeof(char));
    readinfo=(ReadFasta*)malloc((SVM+2)*sizeof(ReadFasta));
    thread=(pthread_t*)malloc(threadnum*sizeof(pthread_t));
    outfile=(FILE **)malloc(threadnum*sizeof(FILE *));
    for(threadno=0; threadno<threadnum; threadno++)
    {
        sprintf(tempstr,"%s/%d.r",workpath,threadno+1);
        outfile[threadno]=fopen(tempstr,"w");
    }
    sprintf(tempstr,"%s/0.fq",workpath);
    fastq=fopen(tempstr,"r");
    //multi process thread
    fileflag=1;
    while(fileflag)
    {
        fileflag=load_fastq(fastq);
        if(readcount%PLL==0)terminalnum=readcount/PLL;
        else terminalnum=readcount/PLL+1;
        if(readcount<=0)break;
        runnumber=0;
        runthreadnum=0;


        pthread_mutex_init(&mutilock,NULL);
        //creat thread
        if(readcount>0)
        {
            for(threadno=0; threadno<threadnum; threadno++)
            {
                threadflag= pthread_create(&thread[threadno], NULL, multithread, NULL);//(multithread)
                if(threadflag)
                {
                    printf("ERROR; return code is %d\n", threadflag);
                    return EXIT_FAILURE;
                }
            }
            //waiting thread
            for(threadno=0; threadno<threadnum; threadno++)pthread_join(thread[threadno],NULL);

        }

        // reference_mapping(1);
    }
    fclose(fastq);
    //clear creat index memory
    free(countin);
    free(databaseindex);
    free(allloc);
    free(REFSEQ);
    free(countin1);
    free(databaseindex1);
    free(allloc1);
    free(read_REFESQ);

    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    fp = fopen("config.txt", "a");
    fprintf(fp, "The Mapping Time: %f sec\n", timeuse);
    fclose(fp);

    for(threadno=0; threadno<threadnum; threadno++)fclose(outfile[threadno]);
    free(outfile);
    free(savework);
    free(readinfo);
    free(thread);
    return 0;
}
