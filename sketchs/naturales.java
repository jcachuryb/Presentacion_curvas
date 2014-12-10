//los puntos fijos que necesitan interpolacion
float[] xD;
float[] yD;
//los componenetes h del polinomio interpolador
float[] hs;
//los aes
float[] as;
//los bes
float[] bs;
//los ces
float[] cs;
//los de
float[] ds;
//----------
void setup(){
   size(600,400);
  xD = new float[4];
  yD = new float[4];
  xD[0]=150;
  yD[0]=300;
  xD[1]=200;
  yD[1]=200;
  xD[2]=350;
  yD[2]=150;
  xD[3]=500;
  yD[3]=200;
  
}

void draw(){
  background(255);
  for(int i=0;i<xD.length;i++){
    fill(255,0,0);
    ellipse(xD[i],yD[i], 5, 5);
  }
  int N=400;
  float delta=(xD[3]-xD[0])/N;
  for(float x=xD[0];x<=xD[3];x=x+delta){
    
    float y=trazado_cubico(xD,yD,x);
    point(x,y);
  
  }
}

float trazado_cubico(float[] xD, float[] yD,float x){
  //declaracion de los coeficientes del 
  //polinomio cubico interpolador
  int n=xD.length;
  hs = new float[n-1];
  as = new float[n];
  bs = new float[n-1];
  cs = new float[n];
  ds = new float[n-1];
  //Matriz de solucion
  float[][] MSol = new float[n][n];
  //vector de solucion para el sistema de ecuaciones
  float[] sVec=new float[n];
  //inicializacion de los coeficientes h
 //a partir de los puntos (xj,yj) 
  for(int i=0;i<xD.length-1;i++){
  
    hs[i]=xD[i+1]-xD[i];
  }
  //inicializacion de los coeficientes h
 //a partir de los puntos yj 
  for(int i=0;i<xD.length;i++){

    as[i]=yD[i];
  }
 //llenar la matriz y vector de soluciones
  MSol[0][0]=1;MSol[n-1][n-1]=1;
  for(int row=0;row<n;row++){
    if(row==0){
      for(int col=1;col<n;col++)
        MSol[row][col]=0;
    }else if(row==n-1){
      for(int col=0;col<n-1;col++)
        MSol[row][col]=0;
    }else{
      for(int col=0;col<n;col++){
        if(abs(row-col)>1){
            MSol[row][col]=0;
        }else if(row==col){
            MSol[row][col]=2*(hs[row-1]+hs[row]);
        }else if(row<col){
            MSol[row][col]=hs[row];
        }else if(row>col){
            MSol[row][col]=hs[row-1];
        }
    
      }
    }
    if(row==0 || row==n-1){
      sVec[row]=0;
    }else{
      sVec[row]=(3/hs[row])*(as[row+1]-as[row])-(3/hs[row-1])*(as[row]-as[row-1]);
    }
    
  } 
  //resolver el sistema de ecuaciones dado SolVec*cs=SVec
  resolver(MSol,sVec,cs,n);
  for(int j=0;j<n-1;j++){
    bs[j]=((as[j+1]-as[j])/hs[j])-(2*cs[j]+cs[j+1])*hs[j]/3;
    ds[j]=(cs[j+1]-cs[j])/(3*hs[j]);
  }
  //buscar el intervalo de x
  int j=0;
  while(x>xD[j+1]){
    j++;
  }
  float aj,bj,cj,dj;
  float hj=x-xD[j];aj=as[j];bj=bs[j];cj=cs[j];dj=ds[j];
  float y=aj + bj*hj + cj*hj*hj+ dj*hj*hj*hj;
return y;

}
void resolver(float a[][], float[] b, float[] x, int n) {
      
  
  
      triangularizar(a, b, n);

      // cálculo de los resultados
      for (int i= n-1; i>=0; i--) {
        float sum= 0.0;
        for (int j= i+1; j<n; j++)
          sum += x[j]*a[i][j];
        x[i]= (b[i]-sum)/a[i][i];
      }
    }    

void triangularizar(float[][] a, float[] b, int n) {
    for(int i=0; i<n; i++) {
      // *
      // dividir fila i por primero de fila i
      float prim= a[i][i];
      for (int j=i; j<n; j++)
        a[i][j] /= prim;
      b[i] /= prim;

      // file j = file j - fila i * primero de fila j
      for(int j= i+1; j<n; j++) {
        prim= a[j][i];
        for(int k= i; k<n; k++)
          a[j][k] -= a[i][k]*prim;
        b[j] -= b[i]*prim;
      }
    }
}