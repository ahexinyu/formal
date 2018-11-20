#include "mecat2ref_defs.h"
#include "output.h"
#include "mecat2ref_aux.h"
#include "../common/diff_gapalign.h"
#include "../common/xdrop_gapalign.h"

#include <algorithm>
using namespace std;

static int MAXC = 0;//MAXC 等于10
static int TECH = TECH_PACBIO;
static int num_output = MAXC;
static const double ddfs_cutoff_pacbio = 0.25;
static const double ddfs_cutoff_nanopore = 0.1;
static double ddfs_cutoff = ddfs_cutoff_pacbio;

static pthread_t *thread;
static int threadnum=2;
static FILE **outfile;
static FILE **refoutfile;
static pthread_mutex_t mutilock; 
static int runnumber=0,runthreadnum=0, readcount,terminalnum;
static int runnumber2=0,runthreadnum2=0,terminalnum2;//*********
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
static int REFcount;//***********
#define FM  200000
float similarity=0;
typedef struct REF_info{
    int refno;
    int reflen;
    char *ref_sequ;
}REF_info;//*************
REF_info *refinfo;//************
char *ref_savework;//************
int mavalue[2000];
int similarity_count;
int read_kmer;
float count_value;
typedef struct {
    char *read_string;
}read_info;
static sim *sc;
static read_info *info;
static int load_ref_f(FILE *fp){
    int readlen,readno,sum=0,flag;
    char *pre;
    REFcount=0;
    pre=ref_savework;
    while((flag=fscanf(fp,"%d\t%d\t%s\n",&readno,&readlen,pre))!=EOF&&REFcount<SVM&&sum<MAXSTR)
    {
        refinfo[REFcount].ref_sequ=pre;
        refinfo[REFcount].refno=readno;
        readlen=strlen(pre);
        refinfo[readcount].reflen=readlen;
        sum=sum+readlen+1;
        pre=pre+readlen+1;
        REFcount++;
    }
    if(flag!=EOF)
        
    {
        refinfo[REFcount].ref_sequ=pre;
        refinfo[REFcount].refno=readno;
        readlen=strlen(pre);
        refinfo[REFcount].reflen=readlen;
        REFcount++;
        return(1);
    }
    else return(0);
}//*
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

static void insert_loc(struct Back_List *spr,int loc,int seedn,float len,long templong)
{//insert_loc(temp_spr,u_k,k+1,BC,sc,templong)
    int list_loc[SI],list_score[SI],list_seed[SI],i,j,minval,mini;int nn=0;int _loc;//在参考基因的位置
    float list_sim[SI];float score_sim[SI];//加了相似度之后的分数
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
    for(i=0;i<SI;i++){score_sim[i]=0;}
    for(i=0;i<SI;i++){
        _loc=(templong*ZV)+list_loc[i];
        nn=(_loc-12)/200;
        printf("n is%d",nn);
        list_sim[i]=(sc[nn].vote);
        score_sim[i]=list_score[i]*list_sim[i];}//考虑相似度
    for(i=0; i<SI; i++)if(minval>score_sim[i])
        {
            minval=score_sim[i];
            mini=i;
        }//找出最低的分数
    if(mini==SM)
    {
        spr->loczhi[SM-1]=loc;
        spr->seedno[SM-1]=seedn;
    }
    else if(mini<SM)
    {
        for(i=mini; i<SM; i++)
        {
            spr->loczhi[i]=list_loc[i+1];
            spr->seedno[i]=list_seed[i+1];
        }
        spr->score--;//删掉最低一个
    }
}
static void insert_loc2(struct Back_List *spr,int loc,int seedn,float len)
{
    int list_loc[SI],list_score[SI],list_seed[SI],i,j,minval,mini,maxi;int maxval;
    for(i=0; i<SM; i++)
    {
        list_loc[i]=spr->loczhi[i];
        list_seed[i]=spr->seedno[i];
        list_score[i]=0;
    }
    list_loc[SM]=loc;
    list_seed[SM]=seedn;
    list_score[SM]=0;
    mini=-1;
    maxi=-1;
    maxval=10;//SM 和SI
    minval=10;
    for(i=0; i<SM; i++)for(j=i+1;j<SI; j++)if(list_seed[j]-list_seed[i]>0&&list_loc[j]-list_loc[i]>0&&fabs((list_loc[j]-list_loc[i])/((list_seed[j]-list_seed[i])*len)-1.0)<ddfs_cutoff)//计算DDF公式
    {
        
        list_score[i]++;
        list_score[j]++;
    }
    for(i=0; i<SI; i++)
    {
        if(maxval<=list_score[i])
        {
            maxval=list_score[i];
            maxi=i;//选出最大的。
        }else if(minval>list_score[i])
        {
            minval=list_score[i];
            mini=i;//选出16里面最小的
        }
    }
    if(maxi==-1)
    {
        for(i=mini; i<SM; i++)
        {
            spr->loczhi[i]=list_loc[i+1];
            spr->seedno[i]=list_seed[i+1];
        }
        spr->score--;
        
    }
    else
    {
        for(i=maxi; i<SM; i++)
        {
            spr->loczhi[i]=list_loc[i+1];
            spr->seedno[i]=list_seed[i+1];
        }
        spr->score--;
    }
}
static void insert_loc3(struct Back_List *spr,int loc,int seedn,float len,long templong)
{//insert_loc(temp_spr,u_k,k+1,BC,templong)
    int list_loc[SI],list_score[SI],list_seed[SI],i,j,minval,mini;int nn=0;int _loc;//在参考基因的位置
    float list_sim[SI];float score_sim[SI];
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
    for(i=0;i<SI;i++){score_sim[i]=0;}
    for(i=0;i<SI;i++){
        _loc=(templong*ZVS)+list_loc[i];
        
        nn=(_loc-12)/200;
        list_sim[i]=(sc[nn].vote);
        score_sim[i]=list_score[i]*list_sim[i];}//考虑相似度
    for(i=0; i<SI; i++)if(minval>score_sim[i])
    {
        minval=score_sim[i];
        mini=i;
    }//找出最低的分数
    if(mini==SM)
    {
        spr->loczhi[SM-1]=loc;
        spr->seedno[SM-1]=seedn;
    }
    else if(mini<SM)
    {
        for(i=mini; i<SM; i++)
        {
            spr->loczhi[i]=list_loc[i+1];
            spr->seedno[i]=list_seed[i+1];
        }
        spr->score--;//删掉最低一个
    }
}
void swapdata(candidate_save *a,candidate_save *b){
    candidate_save temp;
    temp=*a;
    *a=*b;
    *b=temp;
}

void sortdata(candidate_save *can,int num){
    for (int i=0;i<num-1;i++)
    {
        int index=i;
        for(int j=i+1;j<num;j++)
        {
            if(can[j].score>can[i].score)
            {
                index=j;
            }
            if(index!=i)
            {
                swapdata(&can[index], &can[i]);
            }
        }
    }
    
}
static void build_read_index(char *path, char *path1){//buildindex
    unsigned int eit,temp;long start;
    char tempstr[200];
    printf(" path is %s",path);
    sprintf(tempstr, "%s/0.fq",path);
    printf(" tempstr is %s",tempstr);
    int leftnum;
    leftnum=34-2*seed_len;
    FILE *fp;
    fp=fopen(tempstr,"r");
    char *seq;
    int length=get_file_size(path1);//***read文件大小
    printf(" read length is %d",length);
    read_REFESQ=(char *)malloc((length+200000)*sizeof(char));
    seq=read_REFESQ;
    
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
            
        }
        
    }
    
    
    //Max_index
    sumcount1=sumvalue_x(countin1,index_count);
    read_kmer=seqcount1-seed_len+1;
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
            
            
        }
        
    }
    free(read_REFESQ);//****
    free(info);
}

int filter_loc(candidate_save a,candidate_save b){
    int i=0;
    if(fabs((a.loc1-b.loc1)/(a.loc2-b.loc2))>0.7){i=1;return i;}
    else {i=0;return i;}
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
    //sim *sc1;
    //read reference seq
    length=get_file_size(fastafile);
    printf("refernece length is %d",length);
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
    seq[count+1]='\0';
    fclose(fastaindex);
    seqcount=count;
    similarity_count=(seqcount-12)/200+1;
    printf("similarity_count is\n",similarity_count);
    sc=(sim *)malloc((similarity_count+10)*sizeof(sim));
    printf("sim si sucess\n");
    for(int k=0;k<(similarity_count+10);k++){
        sc[k].k_count=0;
        sc[k].simm=0;
        sc[k].TF=0;
        sc[k].LDF=0;//词频
        sc[k].r_count=0;
        sc[k].vote=0;
    }
    //sc1=sc;
    printf("%ld\n",seqcount);
//printf("Constructing look-up table...\n");
    countin=(int *)malloc((indexcount)*sizeof(int));
    for(i=0; i<indexcount; i++)countin[i]=0;
 printf("2222 is suceess\n");
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
    printf("333 is suceess\n");
  int nn=0;//表示区域
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
 printf("444 is suceess\n");
//constructing the look-up table
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
            
            nn=(i-12)/200;//按照200划分，
           if(countin1[eit]>0){
                sc[nn].k_count=sc[nn].k_count+countin1[eit];//在long_read里面出现的次数
            }
            eit=eit<<leftnum;
            eit=eit>>leftnum;
           // printf("%d\n",eit);
        }
    }
    printf("555 is suceess\n");
}

static void get_vote(){
    int eit=0;
    int temp=0;
    int i=0;char *seq;char *readseq;
    readseq=read_REFESQ;//sim *sc1;
    seq=REFSEQ;
    printf("memcy sucuess\n");
    int start=0;//num 有关
    int leftnum=8;int nn=0;
    
    for(int j=0;j<similarity_count-1;j++){
        if(sc[j].k_count>0){
            printf("K-count is %f\n",sc[j].k_count);
            sc[j].LDF=log((read_kmer)/sc[j].k_count);
            printf("LDF is %f\n", sc[j].LDF);
        }
    }
    printf("here is sucuess");
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
            nn=(i-12)/200;
            if(countin[eit]>0){
                //printf("%d\n",cpycount[eit]);
                sc[nn].r_count=sc[nn].r_count+countin[eit];
            } //在参考基因里出现的次数
            eit=eit<<leftnum;
            eit=eit>>leftnum;
        }
    }
     printf("test is sucuess");
    for(int j=0;j<similarity_count;j++){
        if(count_value>0){
           // printf("r_count is %d\n",sc1[j].r_count);
            sc[j].TF=sc[j].r_count/count_value;
        
           // printf("TF is %f\n",sc1[j].TF);}
        
        if((sc[j].TF*sc[j].LDF)>0){
            sc[j].vote=fabs(log(sc[j].TF*sc[j].LDF));}
        else{
            sc[j].vote=1;
        
        };
        //printf("vote is %f\n",sc1[j].vote);
        
    }
    }
        //free(cpycount);
}
int find_location3(int *t_loc,int *t_seedn,int *t_score,long *loc,int k,int *rep_loc,float len,int read_len1, double ddfs_cutoff,long start_loc)//绝了
{
    int i,j,maxval=0,maxi,rep=0,lasti=0;long _loc[200];float list_sim[200];
    for(i=0; i<k; i++){t_score[i]=0;_loc[i]=0;list_sim[i]=0;}
    for(i=0; i<k-1; i++)for(j=i+1; j<k; j++)if(t_seedn[j]-t_seedn[i]>0&&t_loc[j]-t_loc[i]>0&&t_loc[j]-t_loc[i]<read_len1&&fabs((t_loc[j]-t_loc[i])/((t_seedn[j]-t_seedn[i])*len)-1)<ddfs_cutoff)
    {
        t_score[i]++;
        t_score[j]++;
    }
    for(i=0;i<k;i++){_loc[i]=start_loc+t_loc[i];printf(" loc is %d\n",start_loc);printf(" loc is %d\n",_loc[i]);}
    
    int nn=0;
    for(i=0;i<k;i++){
        nn=(_loc[i]-12)/200;
        printf("nn is %d\n",nn);
        list_sim[i]=(sc[nn].vote);
        t_score[i]=t_score[i]*list_sim[i];
    }
    printf("hhhhhhhhhhhhhhhhhhhhhh si sucesss\n");
    for(i=0; i<k; i++)
    {
        if(maxval<t_score[i])
        {
            maxval=t_score[i];
            maxi=i;
            rep=0;
        }
        else if(maxval==t_score[i])
        {
            rep++;
            lasti=i;
        }
    }
    for(i=0; i<4; i++)loc[i]=0;
    if(maxval>=5&&rep==maxval)
    {
        loc[0]=t_loc[maxi],loc[1]=t_seedn[maxi];
        *rep_loc=maxi;
        loc[2]=t_loc[lasti],loc[3]=t_seedn[lasti];
        return(1);
    }
    else if(maxval>=5&&rep!=maxval)
    {
        for(j=0; j<maxi; j++)if(t_seedn[maxi]-t_seedn[j]>0&&t_loc[maxi]-t_loc[j]>0&&t_loc[maxi]-t_loc[j]<read_len1&&fabs((t_loc[maxi]-t_loc[j])/((t_seedn[maxi]-t_seedn[j])*len)-1)<ddfs_cutoff)
        {
            if(loc[0]==0)
            {
                loc[0]=t_loc[j];
                loc[1]=t_seedn[j];
                *rep_loc=j;
            }
            else
            {
                loc[2]=t_loc[j];
                loc[3]=t_seedn[j];
            }
        }
        j=maxi;
        if(loc[0]==0)
        {
            loc[0]=t_loc[j];
            loc[1]=t_seedn[j];
            *rep_loc=j;
        }
        else
        {
            loc[2]=t_loc[j];
            loc[3]=t_seedn[j];
        }
        for(j=maxi+1; j<k; j++)if(t_seedn[j]-t_seedn[maxi]>0&&t_loc[j]-t_loc[maxi]>0&&t_loc[j]-t_loc[maxi]<=read_len1&&fabs((t_loc[j]-t_loc[maxi])/((t_seedn[j]-t_seedn[maxi])*len)-1)<ddfs_cutoff)
        {
            if(loc[0]==0)
            {
                loc[0]=t_loc[j];
                loc[1]=t_seedn[j];
                *rep_loc=j;
            }
            else
            {
                loc[2]=t_loc[j];
                loc[3]=t_seedn[j];
            }
        }
        return(1);
    }
    else return(0);
}


static void reference_mapping(int threadint)
{
    int cleave_num,read_len;
    int mvalue[20000],flag_end;
    long *leadarray,u_k,s_k,loc;
    int count1=0,i,j,k,templong,read_name;//sim *sc1;
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
    //sc1=sc;
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
        if(localnum==terminalnum-1)read_end=readcount;//小于1000条的样子。
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
                                        printf("temp_spr->loczhi[loc-1] is %d\n",u_k);
                                        temp_spr->seedno[loc-1]=k+1;
                                    }
                                    else insert_loc(temp_spr,u_k,k+1,BC,templong);//删除分数最小的。保持在20个左右//
                                    //else insert_loc2(temp_spr,u_k,k+1,BC);
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
                    }//这个时候的score还没有变化
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
                        printf("start_loc is %d",start_loc);
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
                                printf("temp_list[u_k] is %d",temp_list[u_k]);
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
                                //printf("temp_list[u_k] is %d",temp_list[u_k]);
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
                        flag_end=find_location3(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff,start_loc);
                        //flag_end=find_location2(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff);
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
                                if(s_k*1.0/scnt>0.4)temp_spr1->score=0;//有效K_MER
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
                                        //else insert_loc3(temp_spr,u_k,k+1,BC,templong);
                                        else insert_loc2(temp_spr,u_k,k+1,BC);
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
                            //flag_end=find_location(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff,sc,start_loc);
                            flag_end=find_location2(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff);
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
static void reference_map_reference(int threadint)
{
    int  seedlenth=13;int rindexcount=67108864;//***********
    int *table;static long **tableindex1,*tableallloc1,tablesumcount1;//***********
    int leftnum=8;//***********
    int p=0;int pp=0;//**********
    int cleave_num,read_len;
    int mvalue[100000],flag_end;
    long *leadarray,u_k,s_k,loc;
    int count1=0,i,j,k,templong,read_name;
    struct Back_List *database,*temp_spr,*temp_spr1;//
    int repeat_loc = 0,*index_list,*index_spr;
    long location_loc[4],left_length1,right_length1,left_length2,right_length2,loc_list,start_loc;
    short int *index_score,*index_ss;
    int temp_list[200],temp_seedn[200],temp_score[200];
    int localnum,ref_i,read_end,fileid;
    int endnum,ii;
    char *seq,*onedata,onedata1[FM],onedata2[FM],FR;//************
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
    Back_List* rev_database = (Back_List*)malloc(sizeof(Back_List) * j);
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
    u_k = 0;//u_k的数字
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
        localnum=runnumber2;
        runnumber2++;//runnumber是
        pthread_mutex_unlock(&mutilock);
        if(localnum>=terminalnum2)
        {
            fileid=0;
            break;
        }
        if(localnum==terminalnum2-1)read_end=REFcount;//read 条数***********  read_end=readcount
        else read_end=(localnum+1)*PLL;
        for(ref_i=localnum*PLL; ref_i<read_end; ref_i++)
        {
            
            read_name=refinfo[ref_i].refno;
            read_len=refinfo[ref_i].reflen;
            strcpy(onedata1,refinfo[ref_i].ref_sequ);
            
            canidatenum=0;
            for(ii=1; ii<=2; ii++)
            {
                read_len=strlen(onedata1);
                BC=5+(read_len/1000);//BC 是数量
                if(BC>20)BC=20;
                if(ii==1) {
                    onedata=onedata1;
                    index_list = fwd_index_list;
                    index_score = fwd_index_score;
                    database = fwd_database;
                    pnblk = &fnblk;
                }
                /*else if(ii==2){
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
                 }//将序列掉个?
                 /* for(i=0; i<read_len; i++)
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
                 }*/
                endnum=0;
                read_len=strlen(onedata);
                cleave_num=transnum_buchang(onedata,mvalue,&endnum,read_len,seed_len,BC);//切割
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
                        templong=(*leadarray)/ZV;
                        u_k=(*leadarray)%ZV;
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
                                else insert_loc2(temp_spr,u_k,k+1,BC);//insert
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
                for(i=0,index_spr=index_list,index_ss=index_score; i<cc1; i++,index_spr++,index_ss++)if(*index_ss>6)
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
                        if(loc>0)start_loc=(*index_spr-1)*ZVL;
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
                            temp_list[u_k]=temp_spr->loczhi[j]+ZV;
                            temp_seedn[u_k]=temp_spr->seedno[j];
                            u_k++;
                        }
                    }
                    flag_end=find_location2(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff);
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
                    
                    while(canidate_loc[p].score!=0){
                        pp=filter_loc(canidate_loc[p],canidate_temp);
                        if(pp){
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
                                else canidatenum=MAXC;//从高往低排序
                        
                            }
                        else{
                                if(canidate_loc[p].score>=canidate_temp.score){break;}
                                else{
                                    canidate_loc[p]=canidate_temp;
                                    sortdata(canidate_loc,canidatenum);
                                }
                            }
                        p++;
                    }
                    
                    
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
            
            output_results(alns, naln, results, nresults, num_output, refoutfile[threadint]);
            
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
                    } /*else if(ii==2) {
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
                       }*/
                    
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
                                    else insert_loc2(temp_spr,u_k,k+1,BC);
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
                        flag_end=find_location2(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff);
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
                        p=0;
                        pp=0;//insert canidate position or delete this position
                        while(canidate_loc[p].score!=0){
                            pp=filter_loc(canidate_loc[p],canidate_temp);
                            if(pp){
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
                                else canidatenum=MAXC;//从高往低排序
                                
                            }
                            else{
                                if(canidate_loc[p].score>=canidate_temp.score){break;}
                                else{
                                    canidate_loc[p]=canidate_temp;
                                    sortdata(canidate_loc,canidatenum);
                                }
                            }
                            p++;
                        }

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
                
                output_results(alns, naln, results, nresults, num_output, refoutfile[threadint]);
                
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
static void map(char *sonedata,char *sonedata1,TempResult *a,TempResult *b,long length_read,int read_name,FILE *outfile){
    int cleave_num,read_len;
    int mvalue[20000],flag_end;
    long *leadarray,u_k,s_k,loc;
    int count1=0,i,j,k,templong;
    struct Back_List *database,*temp_spr,*temp_spr1;
    int repeat_loc = 0,*index_list,*index_spr;
    long location_loc[4],left_length1,right_length1,left_length2,right_length2,loc_list,start_loc;
    short int *index_score,*index_ss;
    int temp_list[200],temp_seedn[200],temp_score[200];int read_lef,read_right,ref_lef,ref_right;
    int localnum,read_i,read_end,fileid;
    int endnum,ii;int num3,num4;
    char *seq,*onedata,onedata1[RM],onedata2[RM],FR;
    int cc1,canidatenum,loc_seed;
    int num1,num2,BC;
    int low,high,mid,seedcount;
    candidate_save canidate_1,canidate_temp;
    candidate_save  canidate_2;
    seq=REFSEQ;char dir;
    j=seqcount/ZV+5;//z这边考虑下长度要不要改是原来
    
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
    dir=a->read_dir;
    if(dir=='F'){
        strcpy(onedata1,sonedata);
        onedata=onedata1;
    }
    else{
        strcpy(onedata2,sonedata);
        onedata=onedata2;
    }
    canidatenum=0;
    read_len=strlen(onedata);
    BC=5+(read_len/1000);//******这边要改一改。不是1000
    if(BC>20)BC=20;
    if(dir=='F'){
        index_list = fwd_index_list;// int repeat_loc =
        index_score = fwd_index_score;//short int *index_score,*index_ss;
        database = fwd_database;
        pnblk = &fnblk;
    }
   else{
        index_list = rev_index_list;
        index_score = rev_index_score;
        database = rev_database;
        pnblk = &rnblk;
    }
    
    j=0;
    index_list = fwd_index_list;// int repeat_loc =
    index_score = fwd_index_score;//short int *index_score,*index_ss;
    database = fwd_database;
    pnblk = &fnblk;

    endnum=0;
    read_len=strlen(onedata);
    cleave_num=transnum_buchang(onedata,mvalue,&endnum,read_len,seed_len,BC);
    index_spr=index_list;
    index_ss=index_score;
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
                        temp_spr->loczhi[loc-1]=u_k;
                        temp_spr->seedno[loc-1]=k+1;
                    }
                    else insert_loc2(temp_spr,u_k,k+1,BC);
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
   int loc_mid=(a->sb+a->se)/2;
    int block_mid=(loc_mid-12)/ZVL-1;
    for(i=0,index_spr=index_list,index_ss=index_score;i<cc1;i++,index_spr++,index_ss++){
        if(*index_spr==block_mid)
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
            flag_end=find_location2(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len,ddfs_cutoff);
            if(flag_end==0)continue;
            if(temp_score[repeat_loc]<6)continue;//这个改一下待商榷
            canidate_temp.score=temp_score[repeat_loc];
            loc_seed=temp_seedn[repeat_loc];
            location_loc[0]=start_loc+location_loc[0];
            location_loc[1]=(location_loc[1]-1)*BC;//这边要改一下。
            loc_list=location_loc[0];
            ref_lef=location_loc[0]+seed_len-1-a->sb;
            ref_right=a->se-location_loc[0];
            left_length1=location_loc[0]+seed_len-1;
            right_length1=seqcount-location_loc[0];
            read_lef=location_loc[1]+seed_len-1-a->qb;
            read_right=a->qe-location_loc[1];
            left_length2=location_loc[1]+seed_len-1;
            right_length2=read_len-location_loc[1];
            if(ref_lef>=read_lef)num3=read_lef;
            else num3=ref_lef;
            if(ref_right>=read_right)num4=read_right;
            else num4=ref_right;
            if(left_length1>=left_length2)num1=left_length2;
            else num1=left_length1;
            if(right_length1>=right_length2)num2=right_length2;
            else num2=right_length1;
            //缩小范围
            
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
            for(u_k=*index_spr-2,k=num3/ZV,temp_spr1=temp_spr-2; u_k>=0&&k>=0; temp_spr1--,k--,u_k--)if(temp_spr1->score>0)
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
            for(u_k=*index_spr+1,k=num4/ZV,temp_spr1=temp_spr+1; k>0; temp_spr1++,k--,u_k++)if(temp_spr1->score>0)
            {
                start_loc=u_k*ZVL;
                int scnt = min((int)temp_spr1->score, SM);
                for(j=0,s_k=0; j < scnt; j++)if(fabs((start_loc+temp_spr1->loczhi[j]-loc_list)/((temp_spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<ddfs_cutoff)
                {
                    seedcount++;
                    s_k++;
                }
                if(s_k*1.0/scnt>0.4)temp_spr1->score=0;//有效K_MER
            }
            canidate_temp.score=canidate_temp.score+seedcount;//最终分段
            if(dir=='F')canidate_temp.chain='F';
            else canidate_temp.chain='R';
        }
        canidate_1=canidate_temp;
    }
    //这里再来一遍。懂了吗
    naln=0;
    nresults=0;
    int ref_mid=(b->sb+b->se)/2;
    long ref_block_mid=(ref_mid-12)/ZVL-1;
    for(i=0,index_spr=index_list,index_ss=index_score;i<cc1;i++,index_spr++,index_ss++){
        if(*index_spr==ref_block_mid)
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
            flag_end=find_location2(temp_list,temp_seedn,temp_score,location_loc,u_k,&repeat_loc,BC,read_len, ddfs_cutoff);
            if(flag_end==0)continue;
            if(temp_score[repeat_loc]<6)continue;//这个改一下待商榷
            canidate_temp.score=temp_score[repeat_loc];
            loc_seed=temp_seedn[repeat_loc];
            location_loc[0]=start_loc+location_loc[0];
            location_loc[1]=(location_loc[1]-1)*BC;//这边要改一下。
            loc_list=location_loc[0];
            ref_lef=location_loc[0]+seed_len-1-a->sb;
            ref_right=a->se-location_loc[0];
            left_length1=location_loc[0]+seed_len-1;
            right_length1=seqcount-location_loc[0];
            read_lef=location_loc[1]+seed_len-1-a->qb;
            read_right=a->qe-location_loc[1];
            left_length2=location_loc[1]+seed_len-1;
            right_length2=read_len-location_loc[1];
            if(ref_lef>=read_lef)num3=read_lef;
            else num3=ref_lef;
            if(ref_right>=read_right)num4=read_right;
            else num4=ref_right;
            if(left_length1>=left_length2)num1=left_length2;
            else num1=left_length1;
            if(right_length1>=right_length2)num2=right_length2;
            else num2=right_length1;
            
            //find all left seed//
            for(u_k=*index_spr-2,k=num3/ZV,temp_spr1=temp_spr-2; u_k>=0&&k>=0; temp_spr1--,k--,u_k--)if(temp_spr1->score>0)
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
            for(u_k=*index_spr+1,k=num4/ZV,temp_spr1=temp_spr+1; k>0; temp_spr1++,k--,u_k++)if(temp_spr1->score>0)
            {
                start_loc=u_k*ZVL;
                int scnt = min((int)temp_spr1->score, SM);
                for(j=0,s_k=0; j < scnt; j++)if(fabs((start_loc+temp_spr1->loczhi[j]-loc_list)/((temp_spr1->seedno[j]-loc_seed)*BC*1.0)-1.0)<ddfs_cutoff)
                {
                    seedcount++;
                    s_k++;
                }
                if(s_k*1.0/scnt>0.4)temp_spr1->score=0;//有效K_MER
            }
            canidate_temp.score=canidate_temp.score+seedcount;//最终分段
            if(dir=='F')canidate_temp.chain='F';
            else canidate_temp.chain='R';
        }
        canidate_2=canidate_temp;
    }//这里是改之后的
    naln = 0;
    nresults = 0;
    if(canidate_2.score>canidate_1.score){
        extend_candidate(canidate_2,
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
        output_results(alns, naln, results, nresults, num_output, outfile);
        
    }
    else{
        extend_candidate(canidate_1,
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
        output_results(alns, naln, results, nresults, num_output, outfile);//这里的outfile需要改一改
        
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
 int small_meap(TempResult *a,TempResult *b,FILE *upfile){
    
    char onedata[5000],onedata1[5000],onedata2[5000],*onedata3;//onedata2是截取的参考基因部分。
    char raw_read[100000],raw_reference[100000];char *seq;char FR;int j_s=0;
    char f_raw_read[100000];
    int read_id= a->read_id;
    int ref_id=b->read_id;
    int start=a->qb;
    int ref_start=b->qb;
    char dir;
    int read_name;
    read_name=readinfo[read_id].readno;
    if(a->read_dir=='F'){
        strcpy(raw_read,readinfo[read_id].seqloc);
        dir=='F';
    }
    else{
        strcpy(f_raw_read,readinfo[read_id].seqloc);//将read翻转
        onedata3=f_raw_read;
        int len=strlen(f_raw_read);
        int i;
        for(int j=len-1,i=0;j>i;j--,i++){
            FR=f_raw_read[i];
            f_raw_read[i]=f_raw_read[j];
            f_raw_read[j]=FR;
        }
        for(i=0; i<len;i++)
        {
            FR=f_raw_read[i];
            switch(FR)
            {
                case 'A':
                {
                    f_raw_read[i]='T';
                    break;
                }
                case 'T':
                {
                    f_raw_read[i]='A';
                    break;
                }
                case 'C':
                {
                    f_raw_read[i]='G';
                    break;
                }
                case 'G':
                {
                    f_raw_read[i]='C';
                    break;
                }
            }
        }//将字符串掉个
        dir='R';
        printf("faw is %s\n",f_raw_read);
        strcpy(raw_read, f_raw_read);
    }
    strcpy(raw_reference,refinfo[ref_id].ref_sequ);
    /*for(int i=0;i<((a->qe)-(a->qb));i++){
     onedata[i]=raw_read[start+i];
     }
     for(int i=0;i<(b->qe)-(b->qb);i++){
     onedata1[i]=raw_reference[ref_start+i];
     }*/
     //long r_start=(a->sb<b->sb)?a->sb:b->sb;
     //long r_end=(a->se>b->se)?a->se:b->se;
     
     /*for(int i=r_start,j_s=0;i<r_end;i++,j_s++){
     onedata[j_s]=seq[i];
     }*/
     //seq=REFSEQ;//REFSEQ也是原来的，在map函数里面
     //build_small_index(onedata1);//这里应该放read的部分。read和一条参考基因上两个地方比对，考虑了一下，建立两次index吧*/
    long length_read=strlen(raw_read);//这个要加进map函数里面
    map(raw_read,raw_reference,a,b,length_read,read_name,upfile);
    return 0;
    
}



static void* multithread(void* arg)
{
    int localthreadno;int localthreadno2;
    pthread_mutex_lock(&mutilock);
    localthreadno=runthreadnum;
    localthreadno2=runthreadnum2;
    runthreadnum++;
    runthreadnum2++;
    pthread_mutex_unlock(&mutilock);
    reference_mapping(localthreadno);
    reference_map_reference(localthreadno2);
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
	MAXC = maxc;//num_candidate
	TECH = tech;
	num_output = noutput;
    char tempstr[300],fastafile[300];
    char tempstr2[300];//******
    int corenum,readall,refall;//readall是readcount 的条数
    int fileflag,threadno,threadflag;
    FILE *fp,*fastq;
     FILE  *ref_fastq;//*********
    struct timeval tpstart, tpend;
    float timeuse;
    save_work=(char *)malloc((MAXSTR+RM)*sizeof(char));
    fp=fopen("config.txt","r");
    assert(fscanf(fp,"%s\n%s\n%s\n%s\n%s\n%d %d\n%d\n",workpath,fastafile,fastqfile,tempstr,tempstr2,&corenum,&readall,&refall) == 8);//********
     printf("read config is sucess\n");//*********workpath是config的路径
    printf("fastqfile is %s\n",fastqfile);
    fclose(fp);//FASTA是参考基因组文件
    threadnum=corenum;
    
    //building reference index
    gettimeofday(&tpstart, NULL);
    seed_len=13;
    printf("build_read_index pre is sucess\n");
    printf("workpath is %s\n",workpath);
    build_read_index(workpath,fastqfile);//******
    printf("build_read_index sucess\n");
    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    fp = fopen("config.txt", "a");
    fprintf(fp, "The Building read Index Time: %f sec\n", timeuse);
    //fclose(fp);
    //build read_long index
    gettimeofday(&tpstart, NULL);
    printf("get time sis sucess\n");
    creat_ref_index(fastafile);
    printf("get ref_index sucess\n");
    gettimeofday(&tpend, NULL);
    timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
    timeuse /= 1000000;
    //fp = fopen("config.txt", "a");
    fprintf(fp, "The Building  Reference  Index Time: %f sec\n", timeuse);
    fclose(fp);
    get_vote();
    printf("get vote sucess");
    gettimeofday(&tpstart, NULL);

    savework=(char *)malloc((MAXSTR+RM)*sizeof(char));
    ref_savework=(char *)malloc((MAXSTR+RM)*sizeof(char));//********
    readinfo=(ReadFasta*)malloc((SVM+2)*sizeof(ReadFasta));
    refinfo=(REF_info*)malloc((RVM)*sizeof(REF_info));//*********
    thread=(pthread_t*)malloc(threadnum*sizeof(pthread_t));
    outfile=(FILE **)malloc(threadnum*sizeof(FILE *));
    refoutfile=(FILE **)malloc(threadnum*sizeof(FILE *));//*********
    for(threadno=0; threadno<threadnum; threadno++)
    {
        sprintf(tempstr,"%s/%d.r",workpath,threadno+1);
        sprintf(tempstr2,"%s/ref%d.r",workpath,threadno+1);//**********
        
        outfile[threadno]=fopen(tempstr,"w");
        refoutfile[threadno]=fopen(tempstr2,"w");
    }
    sprintf(tempstr,"%s/0.fq",workpath);
    sprintf(tempstr2,"%s/ref.fq",workpath);//******
    fastq=fopen(tempstr,"r");
    ref_fastq=fopen(tempstr2,"r");//***********
    //multi process thread
    fileflag=1;
    while(fileflag)
    {
        fileflag=load_fastq(fastq);
        load_ref_f(ref_fastq);//************8
        if(readcount%PLL==0)terminalnum=readcount/PLL;
        else terminalnum=readcount/PLL+1;
        if(REFcount%PLL==0)terminalnum2=REFcount/PLL;
        else terminalnum2=REFcount/PLL+1;//***********
        if(readcount<=0)break;
        if(REFcount<=0)break;
        runnumber=0;
        runthreadnum=0;
        runnumber2=0;
        runthreadnum2=0;

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
    //free(countin);好像不能free。下面还要用
    free(databaseindex);
    free(allloc);
    //free(REFSEQ);
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

    for(threadno=0; threadno<threadnum; threadno++)
    {
        fclose(outfile[threadno]);
        fclose(refoutfile[threadno]);
    }
    free(outfile);
    free(savework);
    free(save_work);
    free(thread);
    return 0;
}
