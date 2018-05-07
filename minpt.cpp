#define _CRT_SECURE_NO_WARNINGS
#include <random>
#include <algorithm>
#include <numeric>
#include <optional>
#include <fstream>
#include <unordered_map>
#include <filesystem>
#include <queue>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <cstring>
#include <omp.h>
#ifdef _MSC_VER
using namespace std; namespace fs=experimental::filesystem; using I=int;
I bswap(I x){return _byteswap_ulong(x);} string pp(string p){return p;}
#elif defined(__GNUC__)
using namespace std; namespace fs=filesystem; using I=int;
I bswap(I x){return __builtin_bswap32(x);}
string pp(string p){replace(p.begin(),p.end(),'\\','/'); return p;}
#endif
template<class T> using op=optional<T>; using F=double;
using C=char; constexpr F Inf=1e+10,Eps=1e-4,Pi=3.14159265358979323846;
struct V{F x,y,z; V(F v=0):V(v,v,v){} V(F x,F y,F z):x(x),y(y),z(z){}
 F operator[](I i)const{return(&x)[i];} F m()const{return std::max({x,y,z});}};
V operator+(V a,V b){return{a.x+b.x,a.y+b.y,a.z+b.z};}
V operator-(V a,V b){return{a.x-b.x,a.y-b.y,a.z-b.z};}
V operator*(V a,V b){return{a.x*b.x,a.y*b.y,a.z*b.z};}
V operator/(V a,V b){return{a.x/b.x,a.y/b.y,a.z/b.z};}
V operator-(V v){return{-v.x,-v.y,-v.z};} F sq(F v){return v*v;}
V vmin(V a,V b){return{min(a.x,b.x),min(a.y,b.y),min(a.z,b.z)};}
V vmax(V a,V b){return{max(a.x,b.x),max(a.y,b.y),max(a.z,b.z)};}
F dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
V norm(V v){return v/sqrt(dot(v,v));} V refl(V w,V n){return 2*dot(w,n)*n-w;}
V intp(V a,V b,V c,F u,F v) {return a*(1-u-v)+b*u+c*v;}
V cross(V a,V b){return{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
tuple<V,V> odr(V n){F s=copysign(1,n.z),a=-1/(s+n.z),b=n.x*n.y*a;
 V u(1+s*n.x*n.x*a,s*b,-s*n.x),v(b,s+n.y*n.y*a,-n.y); return{u,v};}
op<V> refr(V wi,V n,F et){F t=dot(wi,n),t2=1-et*et*(1-t*t);
 return t2>0?et*(n*t-wi)-n*sqrt(t2):op<V>{};}
struct Rng{mt19937 eng; uniform_real_distribution<F> dist;
 Rng(){}; Rng(I seed){eng.seed(seed); dist.reset();}
 F u(){return dist(eng);} V uD(){F r=sqrt(u()),t=2*Pi*u();
  F x=r*cos(t),y=r*sin(t); return V(x,y,sqrt(max(.0,1-x*x-y*y)));}};
struct Dist{vector<F> c{0}; void add(F v){c.push_back(c.back()+v);}
 void norm(){F sum=c.back(); for(F& v:c){v/=sum;}}
 F p(I i)const{return (i<0||i+1>=I(c.size()))?0:c[i+1]-c[i];}
 I samp(Rng& rn)const{auto it=upper_bound(c.begin(),c.end(),rn.u());
  return clamp(I(distance(c.begin(),it))-1,0,I(c.size())-2);}};
struct Dist2{vector<Dist> ds; Dist m; I w,h;
 void init(const vector<F>& v,I a,I b){w=a; h=b; ds.assign(h,{});
  for(I i=0;i<h;i++){auto& d=ds[i]; for(I j=0;j<w;j++){d.add(v[i*w+j]);}
   m.add(d.c.back()); d.norm();} m.norm();}
 F p(F u,F v)const{I y=min(I(v*h),h-1); return m.p(y)*ds[y].p(I(u*w))*w*h;}
 tuple<F,F> samp(Rng& rn)const{I y=m.samp(rn),x=ds[y].samp(rn);
  return{(x+rn.u())/w,(y+rn.u())/h};}}; struct Ind{I p=-1,t=-1,n=-1;};
namespace T{enum{None=0,AreaL=1<<0,EnvL=1<<1,E=1<<2,D=1<<3,G=1<<4,M=1<<5,
 FresRefl=1<<6,FresTran=1<<7,PRefl=1<<8,Refl=D|G|PRefl|FresRefl,L=AreaL|EnvL,
 Tran=FresTran|M,Fres=FresRefl|FresTran,S=PRefl|Fres|M,NonS=D|G};}
struct Geo{vector<V> ps,ns,ts;}; struct R{V o,d;}; struct Surf{V p,n,t,u,v;
 Surf(){} Surf(V p,V n,V t):p(p),n(n),t(t){tie(u,v)=odr(n);}
 bool op(V wi,V wo)const{return dot(wi,n)*dot(wo,n)<=0;}
 tuple<V,V,V> obn(V wi)const{I i=dot(wi,n)>0;return{i?n:-n,u,i?v:-v};}};
struct H{Surf sp; struct Obj* o;}; struct Tex{I w,h; vector<F> cs,as;
 I fl(I i){I j=i/3,x=j%w,y=j/w; return 3*((h-y-1)*w+x)+i%3;}
 F pf(I i,F e,vector<uint8_t>& ct){return pow(F(ct[fl(i)])/e,2.2);}
 F pf(I i,F e,vector<float>& ct){if(e<0){return ct[fl(i)];}
  auto m=bswap(*(int32_t*)&ct[fl(i)]); return *(float*)&m;}
 template <class T> void loadpxm(vector<F>& c,string p){puts(p.c_str());
  static vector<T> ct; FILE* f=fopen(p.c_str(),"rb"); if(!f){return;}
  F e; fscanf(f,"%*s %d %d %lf%*c",&w,&h,&e); I sz=w*h*3,s=sizeof(T);
  ct.assign(sz,0); c.assign(sz,0); fread(ct.data(),s,sz,f);
  for(I i=0;i<sz;i++){c[i]=pf(i,e,ct);} fclose(f);}
 void loadpfm(string p){loadpxm<float>(cs,p);} void load(string p){
  auto b=fs::path(p); auto pc=b.replace_extension(".ppm").string(); 
  auto pa=(b.parent_path()/fs::path(b.stem().string()+"_alpha.ppm")).string();
  loadpxm<uint8_t>(cs,pc); loadpxm<uint8_t>(as,pa);}
 V ev(V t,bool alpha=0)const{F u=t.x-floor(t.x),v=t.y-floor(t.y);
  I x=clamp(I(u*w),0,w-1),y=clamp(I(v*h),0,h-1),i=w*y+x;
  return alpha?V(as[3*i]):V(cs[3*i],cs[3*i+1],cs[3*i+2]);}};
struct B{V mi=V(Inf),ma=V(-Inf); V center()const{return (mi+ma)*.5;}
 F sa()const{V d=ma-mi; return 2*(d.x*d.y+d.y*d.z+d.z*d.x);}
 V operator[](I i)const{return(&mi)[i];}
 bool isect(const R& r,F tl,F th)const{for(I i=0;i<3;i++){F vd=1/r.d[i];
  F t1=(mi[i]-r.o[i])*vd,t2=(ma[i]-r.o[i])*vd;if(vd<0){swap(t1,t2);}
  tl=max(t1,tl); th=min(t2,th); if(th<tl){return 0;}} return 1;};};
B merge(B b,V p){return{vmin(b.mi,p),vmax(b.ma,p)};}
B merge(B a,B b){return{vmin(a.mi,b.mi),vmax(a.ma,b.ma)};}
F gt(const Surf& s1, const Surf& s2){V d=s2.p-s1.p; F L2=dot(d,d);
 d=d/sqrt(L2); return abs(dot(s1.n,d))*abs(dot(s2.n,-d))/L2;};
struct S{I t; R r; V w; F pcs;}; struct SL{V wo; F d; V fs; F p;};
struct Mat{I t=T::None; V Kd,Ks,Ke; Tex* mapKd=0; F Ni,Ns,an,ax,ay;};
struct Lens{F cr,t,e,ar;}; struct Obj{I t; Mat* M=0; vector<Ind> fs;
 struct{V p,u,v,w; F tf,a,fs,id,ss; vector<Lens> ls; vector<op<B>> pbs;} E;
 struct{Dist dist; F invA;} L; struct{Tex map; F rot; Dist2 dist;} EL;
 op<R> trl(R r)const{F z=0; for(I i=I(E.ls.size())-1;i>=0;i--){
  auto& l=E.ls[i]; z-=l.t; struct LH{F t; V n;}; auto h=[&]()->op<LH>{
   if(l.cr==0){F t=(z-r.o.z)/r.d.z; return t<0?op<LH>{}:LH{t,{}};}
   V c(0,0,z+l.cr),oc=c-r.o; F b=dot(oc,r.d),dt=b*b-dot(oc,oc)+l.cr*l.cr;
   if(dt<0){return{};} F t0=b-sqrt(dt),t1=b+sqrt(dt);
   F t=(r.d.z>0)^(l.cr<0)?min(t0,t1):max(t0,t1); if(t<0){return{};}
   V n=(r.o+t*r.d-c)/l.cr; n=dot(n,-r.d)<0?-n:n; return LH{t,n}; }();
  if(!h){return{};} V p=r.o+r.d*h->t; if(p.x*p.x+p.y*p.y>l.ar*l.ar){return{};}
  r.o=p; if(l.cr==0){continue;} F et=l.e/(i>0&&E.ls[i-1].e!=0?E.ls[i-1].e:1);
  auto wt=refr(-r.d,h->n,et); if(!wt){return{};} r.d=*wt; } return r; }
 F ffd(F id)const{op<R> r; auto& lb=E.ls.back(); for(I i=9;i>0;i--){
  if(r=trl({V(0,0,-lb.t+id),norm(V(lb.ar*i/10,0,-id))})){break;}}
  if(!r){return Inf;} F t=-r->o.x/r->d.x,z=(r->o+t*r->d).z;
  F sz=0; for(auto& l:E.ls){sz+=l.t;} return z<sz?-z:Inf;};
 void initE(string p,V e,V c,V u,F fv,F fd,F fs,F ss,F a){puts(p.c_str());
  E.p=e; E.tf=tan(fv*Pi/180*.5); E.a=a; E.ss = ss; E.fs=fs*.001; 
  E.w=norm(e-c); E.u=norm(cross(u,E.w)); E.v=cross(E.w,E.u); C l[4096];
  ifstream f(p); while(f.getline(l,4096)){if(l[0]=='#'||l[0]=='\0'){continue;}
   F cr,t,eta,ar; sscanf(l,"%lf %lf %lf %lf",&cr,&t,&eta,&ar);
   E.ls.push_back({cr*.001,t*.001,eta,ar*.001*.5});} if(E.ls.empty()){return;}
  F lo=Eps,hi=Inf; for(I i=0;i<99;i++){F mi=(lo+hi)*.5; (ffd(mi)<fd?hi:lo)=mi;}
  E.id=hi; Rng rn(42); I n=64; E.pbs.assign(n,{});
  F cfv=-1,sy=sqrt(E.fs*E.fs/(1+E.a*E.a)),sx=E.a*sy; I m=1<<12;
  for(I i=0;i<n;i++){B b; bool f=0; auto& lb=E.ls.back(); V p(0,0,-lb.t+E.id);
   for(I j=0;j<m;j++){p.x=(i+rn.u())/n*E.fs*.5; F r=sqrt(rn.u()),t=2*Pi*rn.u();
    V pl(r*cos(t)*lb.ar,r*sin(t)*lb.ar,-lb.t); auto rt=trl({p,norm(pl-p)});
    if(rt){f=1; b=merge(b,pl); if(p.x<sx*.5){cfv=max(cfv,rt->d.z);}}}
   if(f){E.pbs[i]=b;}} printf("%lf\n",atan(tan(Pi-acos(cfv))/E.a)*180/Pi*2);}
 void initEL(string p,F rot){EL.map.loadpfm(p); EL.rot=rot*Pi/180;
  auto& cs=EL.map.cs; I w=EL.map.w,h=EL.map.h; vector<F> ls(w*h);
  for(I i=0;i<w*h;i++){V v(cs[3*i],cs[3*i+1],cs[3*i+2]);
   ls[i]=v.m()*sin(Pi*(i/w+.5)/h);} EL.dist.init(ls,w,h);}
 void initAreaL(const Geo& G){for(size_t fi=0;fi<fs.size();fi+=3){
   V a=G.ps[fs[fi].p],b=G.ps[fs[fi+1].p],c=G.ps[fs[fi+2].p],cr=cross(b-a,c-a);
   L.dist.add(sqrt(dot(cr,cr))*.5);} L.invA=1/L.dist.c.back(); L.dist.norm();}
 F D(V wh,V u,V v,V n)const{F x=M->ax,y=M->ay;
  return 1/(Pi*x*y*sq(sq(dot(wh,u)/x)+sq(dot(wh,v)/y)+sq(dot(wh,n))));}
 F G(V wi,V wo,V u,V v,V n)const{auto G1=[&](V w){F c=dot(w,n),s=sqrt(1-c*c);
  F cp=dot(w,u)/s,cs=dot(w,v)/s,a2=sq(cp*M->ax)+sq(cs*M->ay);
  return c==0?0:2/(1+sqrt(1+a2*sq(s/c)));}; return G1(wi)*G1(wo);}
 op<S> samp(Rng& rn,I ty,const Surf& sp,const V& wi)const{
  if(ty&T::E){V rp=2*wi-1; I n=I(E.pbs.size()); auto& lb=E.ls.back();
   if(E.ls.empty()){V d=-norm(V(E.a*E.tf*rp.x,E.tf*rp.y,1));
     return S{T::E,E.p,E.u*d.x+E.v*d.y+E.w*d.z,V(1),1};}
   F sy=sqrt(E.fs*E.fs/(1+E.a*E.a)); V o=rp*V(E.a*sy*.5,sy*.5,0); o.z=E.id;
   F l=sqrt(o.x*o.x+o.y*o.y); I i=clamp(I(l/E.fs*2*n),0,n-1); auto& b=E.pbs[i];
   if(!b){return {};} V bl=b->ma-b->mi; V p=b->mi+bl*V(rn.u(),rn.u(),0);
   F s=l!=0?o.y/l:0,c=l!=0?o.x/l:1; V d=norm(V(c*p.x-s*p.y,s*p.x+c*p.y,p.z)-o);
   auto r=trl({o,d}); if(!r){return{};} F A=bl.x*bl.y,Z=lb.t+E.id;
   F w=d.z*d.z*d.z*d.z*A/(Z*Z); o=r->o; d=r->d; o=E.u*o.x+E.v*o.y+E.w*o.z+E.p;
   d=E.u*d.x+E.v*d.y+E.w*d.z; return S{T::E,o,d,V(w*E.ss),1};}
  else if(ty&T::NonS){auto* mp=M->mapKd; V Kd=(mp?mp->ev(sp.t):M->Kd);
   F wd=Kd.m(),ws=M->Ks.m(); if(wd==0&&ws==0){wd=1;ws=0;}
   F s=wd+ws; wd/=s; ws/=s; ty=rn.u()<wd?T::D:T::G;
   if(ty&T::D){auto[n,u,v]=sp.obn(wi); bool usea=mp&&!mp->as.empty();
    F a=usea?mp->ev(sp.t,1).x:1; V d=rn.uD(); return rn.u()>a ?
     S{T::M,sp.p,-wi,V(1),wd}:S{T::D,sp.p,u*d.x+v*d.y+n*d.z,Kd,wd};}
   else if(ty&T::G){auto[n,u,v]=sp.obn(wi); F u1=rn.u()*2*Pi,u2=rn.u();
    V wh=norm(sqrt(u2/(1-u2))*(M->ax*cos(u1)*u+M->ay*sin(u1)*v)+n);
    V wo=refl(wi,wh); if(sp.op(wi,wo)){return{};}
    return S{T::G,sp.p,wo,ev(T::G,sp,wi,wo)/pdf(T::G,sp,wi,wo),ws};}}
  else if(ty&T::PRefl){return S{T::PRefl,sp.p,refl(wi,sp.n),V(1),1};}
  else if(ty&T::Fres){I i=dot(wi,sp.n)>0;V n=i?sp.n:-sp.n;F et=i?1/M->Ni:M->Ni;
   auto wt=refr(wi,n,et); F Fr=!wt?1:[&](){F cos=i?dot(wi,sp.n):dot(*wt,sp.n);
    F r=(1-M->Ni)/(1+M->Ni); r=r*r; return r+(1-r)*pow(1-cos,5);}();
   return rn.u()<Fr?S{T::FresRefl,sp.p,refl(wi,sp.n),V(1),1}
    : S{T::FresTran,sp.p,*wt,V(et*et),1};} return{};}
 F pdf(I type,const Surf& sp,const V& wi,const V& wo)const{
  if(sp.op(wi,wo)){return 0;} if(type&T::D){return 1/Pi;}
  else if(type&T::G){ V wh=norm(wi+wo); auto[n,u,v]=sp.obn(wi);
   return D(wh,u,v,n)*dot(wh,n)/(4*dot(wo,wh)*dot(wo,n));} return 0;}
 op<SL> sampL(Rng& rn,const Geo& G,const Surf& sp)const{
  if(t&T::AreaL){I i=L.dist.samp(rn); F s=sqrt(max(.0,rn.u()));
   V a=G.ps[fs[3*i].p],b=G.ps[fs[3*i+1].p],c=G.ps[fs[3*i+2].p];
   Surf spL(intp(a,b,c,1-s,rn.u()*s),norm(cross(b-a,c-a)),{});
   V pp=spL.p-sp.p,wo=norm(pp); F p=pdfL(sp,spL,-wo);
   return p==0?op<SL>{}:SL{wo,sqrt(dot(pp,pp)),ev(T::AreaL,spL,{},-wo),p};}
  else if(t&&T::EnvL){auto[u,v]=EL.dist.samp(rn); F t=Pi*v,st=sin(t);
   F p=2*Pi*u+EL.rot; V wo(st*sin(p),cos(t),st*cos(p)); F pL=pdfL(sp,{},-wo);
   return pL==0?op<SL>{}:SL{wo,Inf,ev(T::EnvL,{},{},-wo),pL};} return{};}
 F pdfL(const Surf& sp, const Surf& spL, const V& wo)const{
  if(t&T::AreaL){F G=gt(sp,spL); return G==0?0:L.invA/G;}
  else if(t&&T::EnvL){V d=-wo; F at=atan2(d.x,d.z); at=at<0?at+2*Pi:at;
   F t=(at-EL.rot)*.5/Pi; F u=t-floor(t),v=acos(d.y)/Pi; F st=sqrt(1-d.y*d.y);
   return st==0?0:EL.dist.p(u,v)/(2*Pi*Pi*st*abs(dot(-wo,sp.n)));} return 0;}
 V ev(I ty,const Surf& sp,const V& wi,const V& wo)const{
  if(ty&T::AreaL){return dot(wo,sp.n)<=0 ? V() : M->Ke;}
  else if(ty&T::EnvL){V d=-wo; F at=atan2(d.x,d.z); at=at<0?at+2*Pi:at;
   F t=(at-EL.rot)*.5/Pi; return EL.map.ev({t-floor(t),acos(d.y)/Pi,0});}
  else if(ty&T::G){if(sp.op(wi,wo)){return{};} V wh=norm(wi+wo);
   auto[n,u,v]=sp.obn(wi); V Fr=M->Ks+(1-M->Ks)*pow(1-dot(wo,wh),5);
   return M->Ks*Fr*(D(wh,u,v,n)*G(wi,wo,u,v,n)/(4*dot(wi,n)*dot(wo,n)));}
  else if(ty&T::D){if(sp.op(wi,wo)){return{};}
   auto* mp=M->mapKd; F a=(mp&&!mp->as.empty())?mp->ev(sp.t,1).x:1;
   return (mp?mp->ev(sp.t):M->Kd)*(a/Pi);} return{};}};
struct Scene{Geo G; Obj E; vector<Obj> os; vector<unique_ptr<Tex>> ts;
 vector<Mat> ms; vector<I> Ls; unordered_map<string,I> msmap,tsmap;
 bool useEL=0; Obj* envL(){return useEL?&os.front():nullptr;}
 tuple<const Obj&,F> sampL(Rng& rn)const{I n=I(Ls.size());
  I i=clamp(I(rn.u()*n),0,n-1); return{os.at(Ls[i]),1./n};}
 bool ws(C c){return c==' '||c=='\t';}; F pdfL(){return 1./Ls.size();}
 bool cm(C*& t,const C* c,I n){return !strncmp(t,c,n)&&ws(t[n]);}
 void ss(C*& t){t+=strspn(t," \t");} void sc(C*& t){t+=strcspn(t,"/ \t");}
 F nf(C*& t){ss(t); F v=atof(t); sc(t); return v;}
 V nv(C*& t){V v; v.x=nf(t); v.y=nf(t); v.z=nf(t); return v;}
 I pi(I i,I vn){return i<0?vn+i:i>0?i-1:-1;} Ind pind(C*& t){Ind i; ss(t);
  i.p=pi(atoi(t),I(G.ps.size())); sc(t); if(t++[0]!='/'){return i;}
  i.t=pi(atoi(t),I(G.ts.size())); sc(t); if(t++[0]!='/'){return i;}
  i.n=pi(atoi(t),I(G.ns.size())); sc(t); return i;}
 void nn(C*& t,C name[]){sscanf(t,"%s",name);};
 void loadmtl(string p){ifstream f(p); C l[4096],name[256]; puts(p.c_str());
  while(f.getline(l,4096)){auto* t=l; ss(t); if(cm(t,"newmtl",6)){
    nn(t+=7,name); msmap[name]=I(ms.size()); ms.emplace_back(); continue;}
   if(ms.empty()){continue;} Mat& m=ms.back(); if(cm(t,"Kd",2)){m.Kd=nv(t+=3);}
   else if(cm(t,"Ks",2)) m.Ks=nv(t+=3); else if(cm(t,"Ni",2)) m.Ni=nf(t+=3);
   else if(cm(t,"Ns",2)) m.Ns=nf(t+=3); else if(cm(t,"aniso",5)) m.an=nf(t+=5);
   else if(cm(t,"Ke",2)){m.Ke=nv(t+=3); m.t|=m.Ke.m()>0?T::L:0;}
   else if(cm(t,"illum",5)){ss(t+=6); I v=atoi(t);
    m.t|=v==7?T::Fres:v==5?T::PRefl:T::NonS;}
   else if(cm(t,"map_Kd",6)){nn(t+=7,name); auto it=tsmap.find(name);
    if(it!=tsmap.end()){m.mapKd=ts[it->second].get(); continue;}
    tsmap[name]=I(ts.size()); ts.emplace_back(new Tex);
    ts.back()->load(pp((fs::path(p).remove_filename()/name).string()));
    m.mapKd=ts.back().get();} else{continue;}}}
 void load(string p,string env,F rot){C l[4096],name[256]; ifstream f(p); 
  if(!env.empty()){useEL=1; os.push_back({T::EnvL}); os.back().initEL(env,rot);
  Ls.push_back(0);} while(f.getline(l,4096)){C* t=l; ss(t); I on=I(os.size());
   if     (cm(t,"v" ,1)){G.ps.emplace_back(nv(t+=2));}
   else if(cm(t,"vn",2)){G.ns.emplace_back(nv(t+=3));}
   else if(cm(t,"vt",2)){G.ts.emplace_back(nv(t+=3));}
   else if(cm(t,"f" ,1)){t+=2; if(ms.empty()){ms.push_back({T::D,1});
    os.push_back({T::D,&ms.back()});} Ind is[4]; for(auto& i:is) i=pind(t);
    auto& fs=os.back().fs; fs.insert(fs.end(),{is[0],is[1],is[2]});  
    if(is[3].p!=-1){fs.insert(fs.end(),{is[0],is[2],is[3]});}}
   else if(cm(t,"usemtl",6)){t+=7; nn(t,name); auto& M=ms[msmap.at(name)];
    os.push_back({!M.t?(M.t=T::NonS):M.t,&M}); if(M.t&T::L) Ls.push_back(on);}
   else if(cm(t,"mtllib",6)){nn(t+=7,name);
    loadmtl(pp((fs::path(p).remove_filename()/name).string()));}else continue;}
  for(auto& o:os) if(o.M&&o.M->t&T::AreaL){o.initAreaL(G);}
  for(auto& m:ms){if(!(m.t&T::NonS)){continue;} F r=2/(2+m.Ns);
   F as=sqrt(1-m.an*.9); m.ax=max(1e-3,r/as); m.ay=max(1e-3,r*as);}}
 struct TH{F t,u,v;}; struct Tri{V p1,e1,e2,n; B b; V c; I oi,fi;
  Tri(const V& p1,const V& p2,const V& p3,I oi,I fi)
   : p1(p1),oi(oi),fi(fi){e1=p2-p1; e2=p3-p1; n=norm(cross(p2-p1,p3-p1));
   b=merge(b,p1); b=merge(b,p2); b=merge(b,p3); c=b.center();}
  op<TH> isect(const R& r,F tl,F th)const{V p=cross(r.d,e2),tv=r.o-p1;
   V q=cross(tv,e1); F d=dot(e1,p),ad=abs(d),s=copysign(1,d);
   F u=dot(tv,p)*s,v=dot(r.d,q)*s; if(ad<1e-8||u<0||v<0||u+v>ad) return{};
   F t=dot(e2,q)/d; return t<tl||th<t ? op<TH>{} : TH{t,u/ad,v/ad};}};
 struct N{B b; bool leaf=0; I s,e,c1,c2;}; vector<N> ns;
 vector<Tri> trs; vector<I> ti; void build(){queue<tuple<I,I,I>> q;
  for(size_t oi=0;oi<os.size();oi++){auto& fs=os[oi].fs;
   for(size_t fi=0;fi<fs.size();fi+=3){V a=G.ps[fs[fi].p],b=G.ps[fs[fi+1].p];
    V c=G.ps[fs[fi+2].p]; trs.emplace_back(a,b,c,I(oi),I(fi));}}
  I nt=I(trs.size()); q.push({0,0,nt}); ns.assign(2*nt-1,{}); ti.assign(nt,0);
  iota(ti.begin(),ti.end(),0); mutex mu; condition_variable cv;
  atomic<I> pr=0; I nn=1; bool done=0; auto process=[&](){while(!done){
   auto [ni,s,e]=[&]()->tuple<I,I,I>{unique_lock<mutex> lk(mu);
    if(q.empty()){cv.wait(lk,[&](){return done||!q.empty();});}
    if(done) return{}; auto v=q.front(); q.pop(); return v;}();
   if(done) break; N& n=ns[ni]; for(I i=s;i<e;i++) n.b=merge(n.b,trs[ti[i]].b);
   auto st=[&,s=s,e=e](I ax){auto cmp=[&](I i1, I i2){
    return trs[i1].c[ax]<trs[i2].c[ax];}; sort(&ti[s],&ti[e-1]+1,cmp);};
   auto lf=[&,s=s,e=e](){n.leaf=1; n.s=s; n.e=e; pr+=e-s;if(pr==I(trs.size())){
    done=1;cv.notify_all();}}; if(e-s<2){lf();continue;} F b=Inf; I bi,ba;
   for(I a=0;a<3;a++){thread_local vector<F> l(nt+1),r(nt+1); st(a); B bl,br;
    for(I i=0;i<=e-s;i++){I j=e-s-i; l[i]=bl.sa()*i; r[j]=br.sa()*i;
    bl=i<e-s?merge(bl,trs[ti[s+i]].b):bl;br=j>0?merge(br,trs[ti[s+j-1]].b):br;}
    for(I i=1;i<e-s;i++){F c=1+(l[i]+r[i])/n.b.sa();if(c<b){b=c;bi=i;ba=a;}}}
   if(b>e-s){lf(); continue;} st(ba); I m=s+bi; unique_lock<mutex> lk(mu);
   q.push({n.c1=nn++,s,m}); q.push({n.c2=nn++,m,e}); cv.notify_one();}};
  vector<thread> ths(omp_get_max_threads());
  for(auto& th:ths){th=thread(process);} for(auto& th:ths){th.join();}}
 op<H> isect(const R& r,F tl=Eps,F th=Inf){op<TH> mh,h; I mi,s[99]{},si=0;
  while(si>=0){auto& n=ns.at(s[si--]); if(!n.b.isect(r,tl,th)){continue;}
   if(!n.leaf){s[++si]=n.c1; s[++si]=n.c2; continue;}
   for(I i=n.s;i<n.e;i++) if(h=trs[ti[i]].isect(r,tl,th)){mh=h;th=h->t;mi=i;}}
  if(!mh){return{};} Tri& tr=trs[ti[mi]]; Obj& o=os[tr.oi]; V p=r.o+th*r.d;
  F u=mh->u,v=mh->v; Ind& i=o.fs[tr.fi],j=o.fs[tr.fi+1],k=o.fs[tr.fi+2];
  V n=i.n<0?tr.n:norm(intp(G.ns[i.n],G.ns[j.n],G.ns[k.n],u,v));
  V t=i.t<0?0:intp(G.ts[i.t],G.ts[j.t],G.ts[k.t],u,v); return H{{p,n,t},&o};}};
I main(I,C**v){Scene sc; puts(v[1]); sc.load(v[1],v[2],atof(v[7])); sc.build();
 I ns=atoi(v[5]),mb=atoi(v[6]),w=atoi(v[8]),h=atoi(v[9]),sz=w*h*3;
 vector<V> D(w*h,V()); sc.E.initE(v[3],V(atof(v[10]),atof(v[11]),atof(v[12])),
  V(atof(v[13]),atof(v[14]),atof(v[15])),V(0,1,0),F(atof(v[16])),
  F(atof(v[17])),F(atof(v[18])),F(atof(v[19])),F(w)/h);
 auto L=[&](Rng& rn,V&& uv){V L,T(1),wi=uv; I ty=T::E; Obj* o=&sc.E; Surf sp;
  for(I b=0;b<mb;b++){bool nee=!sc.Ls.empty()&&b>0&&ty&T::NonS;
   auto s=o->samp(rn,ty,sp,wi); if(!s){break;} auto[t,r,w,pcs]=*s; T=T/pcs;
   if(nee){auto[oL,pLs]=sc.sampL(rn); auto sL=oL.sampL(rn,sc.G,sp);
    if(sL){auto[wo,d,fL,pL]=*sL; if(!sc.isect({sp.p,wo},Eps,d*(1-Eps))){
     L=L+T*o->ev(t,sp,wi,wo)*fL/(o->pdf(t,sp,wi,wo)+pL*pLs);}}}
   T=T*w; auto hi=sc.isect(r); if(!hi){hi=H{{},sc.envL()}; if(!hi->o){break;}}
   auto[spH,oH]=*hi; if(oH->t&T::L){ L=L+T*oH->ev(oH->t&T::L,spH,{},-r.d)/
    (b==0||t&T::S?1:oH->pdfL(sp,spH,-r.d)*sc.pdfL()/o->pdf(t,sp,wi,r.d)+1);}
   if(b>3){F q=max(.2,1-T.m()); T=T/(1-q); if(rn.u()<q) break;}
   wi=-r.d; o=oH; sp=spH; ty=o->t&~T::L;} return L;};
 #pragma omp parallel for schedule(dynamic,1)
 for(I i=0;i<w*h;i++){thread_local Rng rn(42+omp_get_thread_num());
  for(I j=0;j<ns;j++) D[i]=D[i]+L(rn,V((i%w+rn.u())/w,(i/w+rn.u())/h,0))/ns;}
 FILE* f=fopen(v[4],"wb"); fprintf(f, "PF\n%d %d\n-1\n",w,h);
 vector<float> d(sz); for(I y=0;y<h;y++) for(I x=0;x<w;x++) for(I i=0;i<3;i++){
  d[3*(y*w+x)+i]=float(D[(h-1-y)*w+(w-1-x)][i]);}
 fwrite(d.data(),4,d.size(),f); fclose(f); return 0;}