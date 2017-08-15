public class testing {
    public static void main(String[] args) {
        int[][] aoi1=new int[5][5];
        aoi1[0]=new int[] {3,7,7,4,6};
        aoi1[1]=new int[] {4,9,7,6,8};
        aoi1[2]=new int[] {8,7,0,7,2};
        aoi1[3]=new int[] {5,0,0,9,5};
        aoi1[4]=new int[] {4,0,6,4,1};
        CSRMatrix m1 = new CSRMatrix(aoi1);
        int[][] aoi2=new int[5][5];
        aoi2[0]=new int[] {4,6,3,6,6};
        aoi2[1]=new int[] {6,0,9,5,1};
        aoi2[2]=new int[] {3,0,7,9,4};
        aoi2[3]=new int[] {0,2,8,4,0};
        aoi2[4]=new int[] {5,0,2,0,8};
        CSRMatrix m2 = new CSRMatrix(aoi2);
        
        print2D(aoi1);
        System.out.println("");
        print2D(aoi2);
        System.out.println("");
        
        CSRMatrix m3 = add(m1,m2);
        print2D(make2DArray(m3));
    }
        
    public static int[] vector_Multiplication(CSRMatrix m, int[] n){
        if(n.length!=m.getNumCols()){
            throw new Error("Wrong dimensions");
        }
        
        int[] output=new int[n.length];
        for(int i=0; i<m.getI().length-1;i++){
            for(int j=m.getI()[i]; j<m.getI()[i+1];j++){
                output[i]+=n[m.getJ()[j]]*m.getD()[j];
            }
        }
        return output;
    }
    
    public static CSRMatrix mult_2(CSRMatrix A, CSRMatrix B){
        int[] Ai=A.getI();
        int[] Aj=A.getJ();
        int[] Ad=A.getD();
        int[] Bi=B.getI();
        int[] Bj=B.getJ();
        int[] Bd=B.getD();
        int[] Ci=new int[A.getNumRows()+1];
        for(int i0=0;i0<Ai.length-1;i0++){
            boolean[] flag=new boolean[B.getNumCols()];
            for(int index1=Ai[i0];index1<Ai[i0+1];index1++){
                int k=Aj[index1];
                for(int index2=Bi[k];index2<Bi[k+1];index2++){
                    int j=Bj[index2];
                    flag[j]=true;
                }
            }
            for(int i=0;i<flag.length;i++){
                if(flag[i]){
                    Ci[i0+1]++;
                }
            }
        }
        for(int i=2;i<Ci.length;i++){
            Ci[i]+=Ci[i-1];
        }
        int[] Cj=new int[Ci[Ci.length-1]];
        int[] Cd=new int[Ci[Ci.length-1]];
        int count=0;
        for(int i0=0;i0<Ai.length-1;i0++){
            int[] temp=new int[B.getNumCols()];
            for(int index1=Ai[i0];index1<Ai[i0+1];index1++){
                int k=Aj[index1];
                for(int index2=Bi[k];index2<Bi[k+1];index2++){
                    int j=Bj[index2];
                    temp[j]+=Ad[index1]*Bd[index2];
                }
            }
            for(int h=0;h<B.getNumCols();h++){
                if(temp[h]!=0){
                    Cj[count]=h;
                    Cd[count]=temp[h];
                    count++;
                }
            }
        }        
        CSRMatrix C = new CSRMatrix(Ci,Cj,Cd);
        return C;
    }
    
    public static CSRMatrix[] decompose(CSRMatrix A){
        int[] diag_di=new int[A.getI().length];
        int[] diag_tempj=new int[diag_di.length-1];
        int[] diag_tempd=new int[diag_di.length-1];
        int diag_count=0;
        
        int[] low_di=new int[A.getI().length];
        int[] low_tempj=new int[(low_di.length*low_di.length)/2];
        int[] low_tempd=new int[(low_di.length*low_di.length)/2];
        int low_count=0;
        
        int[] up_di=new int[A.getI().length];
        int[] up_tempj=new int[(up_di.length*up_di.length)/2];
        int[] up_tempd=new int[(up_di.length*up_di.length)/2];
        int up_count=0;
        
        for(int i=0;i<diag_di.length-1;i++){
            int n0=A.getI()[i];
            int n1=A.getI()[i+1];
            for(int k=n0;k<n1;k++){
                int j=A.getJ()[k];
                if(i>j){
                    low_di[i+1]++;
                    low_tempj[low_count]=j;
                    low_tempd[low_count]=A.getD()[k];
                    low_count++;
                }                
                if(i==j){
                    diag_di[i+1]++;
                    diag_tempj[diag_count]=j;
                    diag_tempd[diag_count]=A.getD()[k];
                    diag_count++;
                }
                if(i<j){
                    up_di[i+1]++;
                    up_tempj[up_count]=j;
                    up_tempd[up_count]=A.getD()[k];
                    up_count++;
                }
            }
        }
        for(int i=0;i<diag_di.length-1;i++){
            diag_di[i+1]+=diag_di[i];
        }
        int[] diag_dd=new int[diag_di[diag_di.length-1]];
        for(int i=0;i<diag_di[diag_di.length-1];i++){
            diag_dd[i]+=diag_tempd[i];
        }
        int[] diag_dj=new int[diag_di[diag_di.length-1]];
        for(int i=0;i<diag_di[diag_di.length-1];i++){
            diag_dj[i]+=diag_tempj[i];
        }CSRMatrix D = new CSRMatrix(diag_di,diag_dj,diag_dd);
        
        
        for(int i=0;i<low_di.length-1;i++){
            low_di[i+1]+=low_di[i];
        }
        int[] low_dd=new int[low_di[low_di.length-1]];
        for(int i=0;i<low_di[low_di.length-1];i++){
            low_dd[i]+=low_tempd[i];
        }
        int[] low_dj=new int[low_di[low_di.length-1]];
        for(int i=0;i<low_di[low_di.length-1];i++){
            low_dj[i]+=low_tempj[i];
        }
        CSRMatrix L = new CSRMatrix(low_di,low_dj,low_dd);
        
        
        for(int i=0;i<up_di.length-1;i++){
            up_di[i+1]+=up_di[i];
        }
        int[] up_dd=new int[up_di[up_di.length-1]];
        for(int i=0;i<up_di[up_di.length-1];i++){
            up_dd[i]+=up_tempd[i];
        }
        int[] up_dj=new int[up_di[up_di.length-1]];
        for(int i=0;i<up_di[up_di.length-1];i++){
            up_dj[i]+=up_tempj[i];
        }
        CSRMatrix U = new CSRMatrix(up_di,up_dj,up_dd);
        
        
        CSRMatrix[] name=new CSRMatrix[3];
        name[0]=D;
        name[1]=L;
        name[2]=U;
        return name;
    }
    
    public static CSRMatrix add(CSRMatrix A, CSRMatrix B){
        int[] ai=A.getI();
        int[] aj=A.getJ();
        int[] ad=A.getD();
        int[] bi=B.getI();
        int[] bj=B.getJ();
        int[] bd=B.getD();
        int count=-1;
        for(int i=0;i<ai.length-1;i++){
            int ak=ai[i];
            int bk=bi[i];
            while(ak<ai[i+1]||bk<bi[i+1]){
                count++;
                if(aj[ak]>bj[bk]){
                    bk++;
                }
                else if(aj[ak]<bj[bk]){
                    ak++;
                }
                else{
                    ak++;
                    bk++;
                }
            }
            while(ak<ai[i+1]){
                count++;
                ak++;
            }
            while(bk<bi[i+1]){
                count++;
                bk++;
            }
        }
        int[] ci=new int[ai.length];
        int[] cj=new int[count];
        int[] cd=new int[count];
        count=-1;
        for(int i=0;i<ai.length-1;i++){
            int ak=ai[i];
            int bk=bi[i];
            while(ak<ai[i+1]&&bk<bi[i+1]){
                count++;
                if(aj[ak]>bj[bk]){
                    ci[i+1]++;
                    cj[count]=bj[bk];
                    cd[count]=bd[bk];
                    bk++;
                }
                else if(aj[ak]<bj[bk]){
                    ci[i+1]++;
                    cj[count]=aj[ak];
                    cd[count]=ad[ak];
                    ak++;
                }
                else{
                    ci[i+1]++;
                    cj[count]=aj[ak];
                    cd[count]=ad[ak]+bd[bk];
                    ak++;
                    bk++;
                }
            }
            while(ak<ai[i+1]){                
                count++;
                ci[i+1]++;
                cj[count]=aj[ak];
                cd[count]=ad[ak];
                ak++;
            }
            while(bk<bi[i+1]){
                count++;
                ci[i+1]++;
                cj[count]=bj[bk];
                cd[count]=bd[bk];
                bk++;
            }
        }
        
        for(int i=0;i<ci.length-1;i++){
            ci[i+1]+=ci[i];
        }
        CSRMatrix C = new CSRMatrix(ci,cj,cd);
        return C;
    }
    
    public static void print2D(int[][] a){
        for(int i=0;i<a.length;i++){
            for(int j=0;j<a[0].length;j++){
                System.out.print(a[i][j]+" ");            
            }
            System.out.println("");
        }
    }
    
    public static int[][] make2DArray(CSRMatrix m){
        int[][] a = new int[m.getNumRows()][m.getNumCols()];
        for(int i=0; i<m.getI().length-1; i++){
            for(int k=m.getI()[i]; k<m.getI()[i+1]; k++){
                a[i][m.getJ()[k]] = m.getD()[k];
            }
        }
        return a;
    }
}
