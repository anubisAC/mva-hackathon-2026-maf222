/* Allele-blind targeted read recruitment, using unique reference-flank k-mers.
 * Input FASTQ is scanned completely. This is not a genome-wide alignment.
 * kseq.h is the MIT-licensed parser from minimap2 v2.31; zlib checks gzip CRCs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#define CAP (1u<<17)
#define BLOOM (1u<<20)
typedef struct {uint64_t key[CAP]; unsigned count[CAP]; unsigned char present[CAP]; unsigned char bloom[BLOOM]; unsigned used;} set;
static set s17,s21;
static int nt(unsigned char c){switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;default:return -1;}}
static unsigned hash(uint64_t x){x^=x>>30;x*=UINT64_C(0xbf58476d1ce4e5b9);x^=x>>27;x*=UINT64_C(0x94d049bb133111eb);x^=x>>31;return (unsigned)x;}
static unsigned bloomkey(uint64_t x){return (unsigned)(x^(x>>20))&(BLOOM-1);}
static unsigned find(set *s,uint64_t key){unsigned i=hash(key)&(CAP-1);while(s->present[i]&&s->key[i]!=key)i=(i+1)&(CAP-1);return i;}
static void insert(set *s,uint64_t key){unsigned i=find(s,key);if(!s->present[i]){if(++s->used>CAP/2){fprintf(stderr,"Too many bait k-mers\n");exit(2);}s->key[i]=key;s->present[i]=1;s->bloom[bloomkey(key)]=1;}}
static void walk(set *s,const char *seq,size_t n,int k,int mode){uint64_t f=0,r=0,mask=(UINT64_C(1)<<(2*k))-1;int valid=0;for(size_t j=0;j<n;j++){int b=nt(seq[j]);if(b<0){f=r=0;valid=0;continue;}f=((f<<2)|b)&mask;r=(r>>2)|((uint64_t)(3-b)<<(2*(k-1)));if(++valid<k)continue;uint64_t key=f<r?f:r;if(mode==0)insert(s,key);else if(s->bloom[bloomkey(key)]){unsigned i=find(s,key);if(s->present[i]&&s->count[i]<2)s->count[i]++;}}}
static int hit(set *s,const char *seq,size_t n,int k){uint64_t f=0,r=0,mask=(UINT64_C(1)<<(2*k))-1;int valid=0;for(size_t j=0;j<n;j++){int b=nt(seq[j]);if(b<0){f=r=0;valid=0;continue;}f=((f<<2)|b)&mask;r=(r>>2)|((uint64_t)(3-b)<<(2*(k-1)));if(++valid<k)continue;uint64_t key=f<r?f:r;if(s->bloom[bloomkey(key)]){unsigned i=find(s,key);if(s->present[i])return 1;}}return 0;}
static void check_gz(gzFile f){int err=0;const char *msg=gzerror(f,&err);if(err!=Z_OK&&err!=Z_STREAM_END){fprintf(stderr,"gzip error: %s\n",msg);exit(2);}}
static int same_name(char *a,char *b){size_t x=strlen(a),y=strlen(b);if(x>2&&a[x-2]=='/'&&(a[x-1]=='1'||a[x-1]=='2'))x-=2;if(y>2&&b[y-2]=='/'&&(b[y-1]=='1'||b[y-1]=='2'))y-=2;return x==y&&memcmp(a,b,x)==0;}
static void emit(FILE *out,kseq_t *r){fprintf(out,"@%s%s%s\n%s\n+\n%s\n",r->name.s,r->comment.l?" ":"",r->comment.l?r->comment.s:"",r->seq.s,r->qual.s);}
static FILE *output(const char *prefix,const char *suffix){char path[4096];snprintf(path,sizeof(path),"%s%s",prefix,suffix);FILE *f=fopen(path,"w");if(!f){perror(path);exit(2);}return f;}
int main(int argc,char **argv){
 if(argc==5&&strcmp(argv[1],"unique")==0){
  gzFile in=gzopen(argv[2],"rb");if(!in)return 2;kseq_t *r=kseq_init(in);int n;
  while((n=kseq_read(r))>=0){walk(&s17,r->seq.s,r->seq.l,17,0);walk(&s21,r->seq.s,r->seq.l,21,0);}check_gz(in);if(n!=-1)return 2;kseq_destroy(r);gzclose(in);
  in=gzopen(argv[3],"rb");if(!in)return 2;r=kseq_init(in);
  while((n=kseq_read(r))>=0){walk(&s17,r->seq.s,r->seq.l,17,1);walk(&s21,r->seq.s,r->seq.l,21,1);}check_gz(in);if(n!=-1)return 2;kseq_destroy(r);gzclose(in);
  FILE *out=fopen(argv[4],"w");if(!out)return 2;unsigned count17=0,count21=0;
  for(unsigned i=0;i<CAP;i++){if(s17.present[i]&&s17.count[i]==1){fprintf(out,"17\t%llu\n",(unsigned long long)s17.key[i]);count17++;}if(s21.present[i]&&s21.count[i]==1){fprintf(out,"21\t%llu\n",(unsigned long long)s21.key[i]);count21++;}}
  fclose(out);printf("Unique k-mers: k17=%u/%u, k21=%u/%u\n",count17,s17.used,count21,s21.used);return count17&&count21?0:2;
 }
 if(argc!=6&&argc!=7){fprintf(stderr,"Usage: recruit_reads unique masked_baits.fa genome.fa.gz kmers.tsv\n       recruit_reads scan kmers.tsv r1.fq.gz r2.fq.gz output_prefix [max_pairs]\n");return 2;}
 if(strcmp(argv[1],"scan"))return 2;FILE *input=fopen(argv[2],"r");if(!input)return 2;int k;unsigned long long key;
 while(fscanf(input,"%d\t%llu",&k,&key)==2){if(k!=17&&k!=21)return 2;insert(k==17?&s17:&s21,key);}if(!feof(input))return 2;fclose(input);
 gzFile a=gzopen(argv[3],"rb"),b=gzopen(argv[4],"rb");if(!a||!b)return 2;gzbuffer(a,1048576);gzbuffer(b,1048576);kseq_t *one=kseq_init(a),*two=kseq_init(b);
 FILE *o1=output(argv[5],"_R1.fastq"),*o2=output(argv[5],"_R2.fastq"),*names=output(argv[5],"_recruitment.tsv");fprintf(names,"read_name\tk17_r1\tk17_r2\tk21_r1\tk21_r2\n");
 uint64_t total=0,kept=0,strict=0,limit=argc==7?strtoull(argv[6],NULL,10):0;int x,y;double start=(double)time(NULL);int complete=0;
 while(1){x=kseq_read(one);y=kseq_read(two);if(x==-1&&y==-1){complete=1;break;}if(x<0||y<0||one->qual.l!=one->seq.l||two->qual.l!=two->seq.l||!same_name(one->name.s,two->name.s)){fprintf(stderr,"Malformed or mismatched FASTQ at pair %llu (%d,%d)\n",(unsigned long long)total+1,x,y);return 2;}total++;
  int h1=hit(&s17,one->seq.s,one->seq.l,17),h2=hit(&s17,two->seq.s,two->seq.l,17);
  if(h1||h2){int q1=hit(&s21,one->seq.s,one->seq.l,21),q2=hit(&s21,two->seq.s,two->seq.l,21);emit(o1,one);emit(o2,two);fprintf(names,"%s\t%d\t%d\t%d\t%d\n",one->name.s,h1,h2,q1,q2);kept++;if(q1||q2)strict++;}
  if(total%5000000==0){fprintf(stderr,"pairs=%llu recruited=%llu compressed_read_bytes=%lld elapsed_s=%.0f\n",(unsigned long long)total,(unsigned long long)kept,(long long)(gzoffset(a)+gzoffset(b)),(double)time(NULL)-start);fflush(stderr);}
  if(limit&&total>=limit)break;
 }
 check_gz(a);check_gz(b);if(fclose(o1)||fclose(o2)||fclose(names))return 2;
 FILE *summary=output(argv[5],"_summary.json");fprintf(summary,"{\"pairs_scanned\":%llu,\"pairs_recruited_k17\":%llu,\"pairs_recruited_k21\":%llu,\"complete_scan\":%s,\"elapsed_seconds\":%.0f,\"gzip_streams_valid_to_end\":%s}\n",(unsigned long long)total,(unsigned long long)kept,(unsigned long long)strict,complete?"true":"false",(double)time(NULL)-start,complete?"true":"false");fclose(summary);kseq_destroy(one);kseq_destroy(two);gzclose(a);gzclose(b);return 0;
}
