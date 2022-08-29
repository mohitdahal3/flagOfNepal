//import java.lang.Math;
//import java.util.Arrays;


float s = 400;
PVector A;
PVector B;
PVector C;
PVector D;
PVector E;
PVector F;
PVector G;
PVector H;
PVector I;
PVector J;
PVector K;
PVector L;
PVector M;
PVector N;
PVector O;
PVector i1;
PVector P;
PVector Q;
PVector R;
PVector S;
PVector T;
PVector U;
PVector V;
PVector W;

void setup(){
  size(600 , 500);
  float temp;
  A = new PVector(0,0);
  B = new PVector(0,s);
  temp = (s*3)/(4);
  C = new PVector(temp , s);
  D = new PVector(0 , s - (dis(B,C)));
  E = vec(C,D);
  E.setMag(dis(B,C));
  E = PVector.add(C , E);
  F = new PVector(0 , E.y);
  G = new PVector(dis(B,C) , E.y);
  H = new PVector( (dis(B,C)) / (4) , B.y);
  I = lineLineIntersection(A,G,H,new PVector(H.x , 0));
  J = new PVector(A.x , (dis(A , F))/(2) );
  K = lineLineIntersection(A ,G , J , new PVector(s , J.y));
  L = lineLineIntersection(J,K,H,I);
  M = lineLineIntersection(J,G,H,I);
  N = vec(M , H);
  N.setMag(disPointToLine(D,C,M));
  N = PVector.add(M , N);
  O = new PVector(0 , M.y);
  i1 = lineLineIntersection(O,M,A,G);
  P = vec(M,O);
  P.setMag(sqrt( (dis(L,N) * dis(L,N)) - (dis(L,M) * dis(L,M)) ));
  P = PVector.add(M , P);
  
  Q = vec(O,M);
  Q.setMag(sqrt( (dis(L,N) * dis(L,N)) - (dis(L,M) * dis(L,M)) ));
  Q = PVector.add(M , Q);
  R = circleCirlceIntersection(L,dis(L,N) , dis(M,N))[1];
  S = circleCirlceIntersection(L , dis(L,N) , dis(M,N))[0];
  T = lineLineIntersection(R , S , H , I);
  U = new PVector(0,B.y - (dis(B,F)) / (2));
  V = lineLineIntersection(C,E,U , new PVector(U.x + 100 , U.y));
  W = lineLineIntersection(U,V,H,I);
  
}

void draw(){
  background(200);
  translate(140 , 50);
  stroke(0);
  strokeWeight(2);
  noFill();
  dline(A,B);
  dline(B,C);
  dline(C,D);
  dline(F,G);
  dline(A,G);
  
  beginShape();
    vertex(A.x,A.y);
    vertex(B.x , B.y);
    vertex(C.x , C.y);
    vertex(E.x , E.y);
    vertex(G.x , G.y);
  endShape(CLOSE);
  
  dline(H,I);
  dline(J,K);
  dline(G,J);
  dline(O,i1);
  float theta;
  theta = asin( (dis(L,M))/(dis(L,N)) );
  arc(L.x , L.y , 2*dis(L,N) , 2*dis(L,N) , theta , PI - theta);
  arc(M.x , M.y , 2*dis(M,Q) , 2*dis(M,Q) , 0 , PI);
  arc(N.x , N.y , 2*dis(M,N) , 2*dis(M,N) , PI + PI/8 , TWO_PI - PI/8);
  arc(T.x , T.y , 2*dis(T,S) , 2*dis(T,S) , PI , TWO_PI);
  
  dline(R,S);
  dline(U,V);
  circle(W.x,W.y,2*dis(M,N));
  circle(W.x,W.y,2*dis(L,N));
  
  drawMoonSpikes();
  drawSunSpikes();
  noLoop();
}

void dline(PVector a , PVector b){
  line(a.x , a.y , b.x , b.y);
}

float dis(PVector a , PVector b){
  return a.copy().sub(b).mag(); 
}


PVector vec(PVector a , PVector b){
  return b.copy().sub(a);
}

void poi(PVector a){
  strokeWeight(6);
  point(a.x , a.y);
}


PVector lineLineIntersection(PVector A, PVector B, PVector C, PVector D){
  // Line AB represented as a1x + b1y = c1
  float a1 = B.y - A.y;
  float b1 = A.x - B.x;
  float c1 = a1*(A.x) + b1*(A.y);
      
  // Line CD represented as a2x + b2y = c2
  float a2 = D.y - C.y;
  float b2 = C.x - D.x;
  float c2 = a2*(C.x)+ b2*(C.y);
      
  float determinant = a1*b2 - a2*b1;
  
  
  if (determinant == 0){
    
    // The lines are parallel. This is simplified
    // by returning a pair of FLT_MAX
    return new PVector(0,0);
  }else{
    float x = (b2*c1 - b1*c2)/determinant;
    float y = (a1*c2 - a2*c1)/determinant;
    return new PVector(x, y);
  }
}


float disPointToLine(PVector a , PVector b , PVector c){
  float x1 = a.x;
  float y1 = a.y;
  float x2 = b.x;
  float y2 = b.y;
  float x3 = c.x;
  float y3 = c.y;
  float m = (y2-y1)/(x2-x1);
  float num = (m*x3) - (y3) - (m*x1) + (y1);
  float den = sqrt((m*m) + 1);
  return abs(num/den);
}


PVector[] circleCirlceIntersection(PVector c1 , float r1 , float r2){
  float y = ( (2 * r1 * r1) - (r2 * r2) ) / (2*r1);
  float x1 = sqrt( (r1*r1) - (y*y));
  float x2 = x1 * -1;
  y += c1.y;
  x1 += c1.x;
  x2 += c1.x;
  PVector[] arr = {new PVector(x1,y) , new PVector(x2,y)};
  return arr;
}


void drawMoonSpikes(){
  pushMatrix();
  translate(N.x , N.y);
  float r = dis(M , N);
  float a = asin((N.y-S.y)/(r));
  float theta = (PI - 2*a) / 8;
  float t = -a;
  noFill();
  strokeWeight(2);
  beginShape();
   for(int i = 0; i < 9; i++){
     float x = r * cos(t);
     float y = r * sin(t);
     float x1 = r * cos(t-theta);
     float y1 = r * sin(t-theta);
     t -= theta;
     PVector tvec = new PVector(0 , -1 * dis(N,T));
     PVector mid = new PVector((x+x1)/(2) , (y+y1)/(2));
     PVector far = vec(tvec , mid);
     far.setMag(dis(T,S));
     far = PVector.add(far , tvec);
     vertex(x,y);
     if(i!=8){
       vertex(far.x , far.y);
     }
     
   }
   endShape();
   popMatrix();
}


void drawSunSpikes(){
  pushMatrix();
  translate(W.x , W.y);
  float t = 0;
  float r = dis(M,N);
  float theta = TWO_PI / 12;
  noFill();
  strokeWeight(2);
  beginShape();
  for (int i = 0; i < 12; i++){
    float x = r * cos(t);
    float y = r * sin(t);
    float x1 = r * cos(t + theta);
    float y1 = r * sin(t + theta);
    PVector mid = new PVector((x+x1)/(2) , (y+y1)/(2));
    PVector far = vec(new PVector(0,0) , mid);
    far.setMag(dis(L,N));
    t += theta;
    vertex(x,y);
    vertex(far.x , far.y);
  }
  endShape(CLOSE);
  
  popMatrix();
}
