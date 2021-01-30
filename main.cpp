//
//  main.cpp
//  CPPcode
//
//  Created by 王建凯 on 2021/1/27.
//  Copyright © 2021 王建凯. All rights reserved.
//

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<cstdlib>
#include<time.h>  //clock()
using namespace std;
//创建动态数组
int *aarray1(int m)
{
    int *p;
    p=(int *)malloc((m)*sizeof(int));
    return p;
}
//释放二维数组所占用的内存
void freearray1(int *p)
{
    free(p);
}
float **aarray(int m,int n)
{
    int i;
    float **p;
    p=(float **)malloc((m)*sizeof(float *));
    p[0]=(float *)malloc((m)*(n)*sizeof(float));
    for(i=1;i<(m);i++)    p[i]=p[i-1]+(n);
    return p;
}
//释放二维数组所占用的内存
void freearray(float **p)
{
    free(*p);
    free(p);
}
//用来计算两个向量坐标的累加和
float sumofsquare(float **aaa,float **dataSetArray,int columns,int ii,int jj){
//    cout<<"aaa "<<aaa[0][0]<<"  "<<aaa[0][1]<<endl;
//    cout<<"aaa "<<aaa[1][0]<<"  "<<aaa[1][1]<<endl;
//    cout<<"aaa "<<aaa[2][0]<<"  "<<aaa[2][1]<<endl;
//    cout<<"aaa "<<aaa[3][0]<<"  "<<aaa[3][1]<<endl;
    int k;
    float ss1;
    float sum333 = 0;
    for(k=0;k<columns;k++){
        ss1 = (aaa[ii][k]-dataSetArray[jj][k])*(aaa[ii][k]-dataSetArray[jj][k]);
        //                    cout<<ss1<<endl;
//        cout<<"计算平方和ss1 "<<ss1<<endl;
//        cout<<"计算方差的坐标 "<<dataSetArray[jj][k]<<endl;
        sum333 += ss1;
    }
    
    return sum333;

}
//聚类后平方误差计算
float gate(int *means,float **aaa,float **dataSetArray,int value_of_k,int rows,int columns){
    int i=0,j=0;
    float sum = 0;
   
//    float a[value_of_k];
//    float sum = 0,sum1 = 0;
    float sum22 = 0;
    while(i<value_of_k){
        float sum11 = 0;
//        for(j=0;j<rows;j++){
        while(j<rows){
                if(means[j] == i){
                sum = sumofsquare(aaa,dataSetArray,columns,i,j);
                sum11 += sum;
//                cout<<"sum11 "<<sum11<<endl;
                j=j+1;
                }else{
                    j= j+1;
                }
        }
        sum22 = sum22+sum11;
//        cout<<"sum22 "<<sum22<<endl;
        i=i+1;
     }
    return sum22;
}
float **getmean(int *means,float **dataSetArray,int value_of_k,int rows,int columns){//得到四个簇的聚类中心
    float sum= 0;
//    float sumX=0;
    int i = 0;
    //    int k=0,kk=0,kkk=0,kkkk=0;
    //    for(int k=0;k<*value_of_k;k++){
    float **p1 = aarray(value_of_k,columns);
    int kk1 = 0;
    int *a = aarray1(columns);
    //        float sumXX=0.0,sumYY=0.0;
    //        float sumXXX=0.0,sumYYY=0.0;
    //        float sumXXXX=0.0,sumYYYY=0.0;
    //        float **p1 = new float[*rows][*columns];
    for(int i3=0;i3<value_of_k;i3++){
        for(int i2=0;i2<columns;i2++){
            sum = 0;
            kk1 = 0;
            i = 0;
            while(i<rows){
                if(means[i] == i3){
                    kk1=kk1+1;
                    sum += dataSetArray[i][i2];
                    i = i+1;
                }else{
                    i = i+1;
                    continue;
                }
                
            }
            if(kk1 != 0){
                a[i2]=sum/kk1;
            }
        
        }
        for(int i4=0;i4<columns;i4++){
            p1[i3][i4]=a[i4];
            //cout<<"打印P数组"<<p1[i3][i4]<<endl;
        }
    }
    return p1;
}

float getdistance(int i,int j,float **dataSetArray,int value_of_k,int rows,int columns){//得到数据集数组中第i个点和第j个点的欧式距离
    //    loadDataSet();
//    float d;
    float sum = 0;
    float dd;
    //    d = sqrt((dataSetArray[i][0]-dataSetArray[j][0])*(dataSetArray[i][0]-dataSetArray[j][0]) + (dataSetArray[i][1]-dataSetArray[j][1])*(dataSetArray[i][1]-dataSetArray[j][1]));
    for(int i1= 0;i1<columns;i1++){
        dd = (dataSetArray[i][i1]-dataSetArray[j][i1])*(dataSetArray[i][i1]-dataSetArray[j][i1]);
        sum+=dd;
    }
    return sqrt(sum);
}
int minmal(int i,float **distance,int value_of_k,int rows,int columns){//为每个点计算最近的聚类中心，distance[80][4]中共有80个点，每个点距离4个聚类中心的4个距离 求最小
    int n = 0;
    float min = distance[i][0];
    int k =0;
    while(n<value_of_k){
        if(min <= distance[i][n]){
            n=n+1;
        }else{
            min = distance[i][n];
            k=n;
            n=n+1;
        }
        
    }
    return k;//返回最小下标
}
float getdistance1(int i,int j,float **dataSetArray,float **p,int value_of_k,int rows,int columns){//求 dataSetArray[80][2]第i个点与 float p[4][2]中第j 个点的距离
//    float d;
    //    d = sqrt((dataSetArray[i][0]-p[j][0])*(dataSetArray[i][0]-p[j][0]) + (dataSetArray[i][1]-p[j][1])*(dataSetArray[i][1]-p[j][1]));
    float sum = 0;
    float dd;
    for(int i1= 0;i1<columns;i1++){
        dd = (p[j][i1]-dataSetArray[i][i1])*(p[j][i1]-dataSetArray[i][i1]);
        sum+=dd;
    }
    return sqrt(sum);
}
int main(){
    FILE *fpp;
    FILE *fpp1;
    int oo =1;
    int kvalue=0;
    int numofcols = 0;
    int numofrows = 0;
    int *value_of_k = &kvalue;
    int *rows = &numofrows;
    int *columns = &numofcols;
//    int *value_of_k,*rows,*columns;//三个参数
    if((fpp=fopen("/Users/wangjiankai/Documents/CPPcode/CPPcode/testSet1.txt","r"))==NULL)    fprintf(stderr,"cannot open data.txt!\n");
    if(fscanf(fpp,"%d %d %d",value_of_k,rows,columns)!=3)        fprintf(stderr,"load error!\n");
    //    float cents[4][2];
    //    double (*pp)[2];
    clock_t start,end;
    cout<<"kvalue"<<kvalue<<endl;
    float **ppp =  aarray(kvalue,numofcols);
    string str,line,x;
    float **dataSetArray = aarray(numofrows,numofcols);//从文件中读取数据到这个数组
    float a,b,c;
    float **distance = aarray(numofrows,kvalue);//记录80个点到4个聚类中心的距离
    int j=0,m,k,f;
    int *means = aarray1(numofrows);//记录每个点到最近的聚类中心
    float sumX=0.0,sumY=0.0;
    //    float cent[4][2];
//    float xx,yy;
//    int kk=0;
    //    FILE *fout;
    ofstream fout("/Users/wangjiankai/Documents/CPPcode/CPPcode/testSet11.txt");
    //    ifstream fin("D:\\c++程序\\testSet2.txt");
    //    fin.open("D://testSet",ios::out);
    //    if(fin==NULL){
    //        cout<<"文件不能打开"<<endl;
    //    }
    //    else{
    //        cout<<"文件打开"<<endl;
    //    }
    if((fpp1=fopen("/Users/wangjiankai/Documents/CPPcode/CPPcode/testSet2.txt","r"))==NULL)    fprintf(stderr,"cannot open data.txt!\n");
    for(int i=0;i<numofrows;i++)
        for(int j=0;j<numofcols;j++)
            fscanf(fpp1,"%f",&dataSetArray[i][j]);  //读取数据点
    
    
    //    while(getline(fin,line)){
    //        istringstream stream(line);
    //        stream>>a>>b;
    ////        cout<<a<<" "<<b<<endl;
    //        dataSetArray[i][0]=a;
    //        dataSetArray[i][1]=b;
    ////        fprintf(fout,"%f %f\n",dataSetArray[i][0],dataSetArray[i][1]);
    //        i=i+1;
    //        j=j+1;
    //    }
    start = clock();
    //    c=getdistance(1,2,dataSetArray);
    //    cout<<c<<endl;
    for(int k1=0;k1<numofrows;k1++){
        for(m=0;m<kvalue;m++){
            c = getdistance(k1,m,dataSetArray,kvalue,numofrows,numofcols);
//        cout<<c<<endl;
            distance[k1][m] = c;
        }
    cout<<"点与对应中心的距离"<<k1<<"\t"<<distance[k1][0]<<"\t"<<distance[k1][1]<<"\t"<<distance[k1][2]<<"\t"<<distance[k1][3]<<endl;// 输出距离
    }
    
    for(int l =0;l<numofrows;l++){
        if(l < kvalue) {
            means[l]=l;//数据集中选择前k个点作为初始聚类中心
        }else{
            f = minmal(l,distance,kvalue,numofrows,numofcols);
            means[l]=f;
        }
    }
    for(int i=1;i<=numofrows;i++){
        cout<<"第"<<i-1<<"个点，对应中心为"<<means[i-1]<<endl;
    }
    //重新计算质心
    
    float **pp = aarray(kvalue,numofcols);                // 第一次计算质心 ，选取前4个点
    for(int i=0;i<kvalue;i++){
        for(int j=0;j<numofcols;j++){
            pp[i][j]=dataSetArray[i][j];
        }
    }
    cout<<"pp"<<pp[0][0]<<"  "<<pp[0][1]<<endl;
    cout<<pp[1][0]<<"  "<<pp[1][1]<<endl;
    cout<<pp[2][0]<<"  "<<pp[2][1]<<endl;
    cout<<pp[3][0]<<"  "<<pp[3][1]<<endl;
    float sum = gate(means,pp,dataSetArray,kvalue,numofrows,numofcols);//第一次平方差计算
    cout<<"第一次方差计算："<<sum<<endl;
    ppp = getmean(means,dataSetArray,kvalue,numofrows,numofcols);//将坐标求平均，重新计算聚类中心
    cout<<"ppp"<<ppp[0][0]<<"  "<<ppp[0][1]<<endl;
    cout<<ppp[1][0]<<"  "<<ppp[1][1]<<endl;
    cout<<ppp[2][0]<<"  "<<ppp[2][1]<<endl;
    cout<<ppp[3][0]<<"  "<<ppp[3][1]<<endl;
    for(int k2=0;k2<numofrows;k2++){
        for(m=0;m<kvalue;m++){
            c = getdistance1(k2,m,dataSetArray,ppp,kvalue,numofrows,numofcols);
//                       cout<<"C值"<<" "<<c<<endl;
            distance[k2][m] = c;//重新得到距离数组
        }
                cout<<"重新得到距离数组"<<distance[k2][0]<<" "<<distance[k2][1]<<" "<<distance[k2][2]<<" "<<distance[k2][3]<<endl;//输出
    }
    for(int l =0;l<numofrows;l++){
        f = minmal(l,distance,kvalue,numofrows,numofcols);//得到每个点的聚类中心，存储在数组means
        cout<<"聚类中心means"<<l<<" "<<f<<endl;
        means[l]=f;
        
    }
    float sum1 = gate(means,ppp,dataSetArray,kvalue,numofrows,numofcols); //第二次计算的最小平方误差
    //    double xxx = fabs(sum1-sum);
    //    cout<<xxx<<endl;
    //    ppp = getMean(means,dataSetArray);
    cout<<sum1<<endl;
    while(fabs(sum-sum1)>0.5){
        sum = sum1;
        cout<<"前后误差的差"<<fabs(sum-sum1)<<endl;
        ppp = getmean(means,dataSetArray,kvalue,numofrows,numofcols);//循环计算聚类中心 存储在ppp[4][2]中
        cout<<ppp[0][0]<<"  "<<ppp[0][1]<<endl;
        cout<<ppp[1][0]<<"  "<<ppp[1][1]<<endl;
        cout<<ppp[2][0]<<"  "<<ppp[2][1]<<endl;
        cout<<ppp[3][0]<<"  "<<ppp[3][1]<<endl;
        for(int k3=0;k3<numofrows;k3++){
            for(m=0;m<kvalue;m++){
                c = getdistance1(k3,m,dataSetArray,ppp,kvalue,numofrows,numofcols);
//                              cout<<c<<endl;
                distance[k3][m] = c;//更新距离数组 distance
            }
        }
        for(int l =0;l<numofrows;l++){
            f = minmal(l,distance,kvalue,numofrows,numofcols);
            cout<<l<<" "<<f<<endl;
            means[l]=f;    //每个点所属聚类中心下标存储在means中
        }
        
        float sum2 = gate(means,ppp,dataSetArray,kvalue,numofrows,numofcols);//计算聚类误差平方和
        //        ppp = getMean(means,dataSetArray);//重新计算质心
        //        for(int ij= 0;ij<*value_of_k;ij++){
        //        cout<<ij<<"   "<<ppp[ij][0]<<"  "<<ppp[ij][1]<<endl;
        //        }
        sum1=sum2;
        
    }
    
    for(int i= 0;i < numofrows;i++){
        for(int i5 =0;i5<numofcols;i5++){
            fout<<dataSetArray[i][i5]<<" ";
            if(i5 == numofcols -1)     fout<<means[i]<<endl;
        }
        
    }
    //    fin.close();
    fout.close();
    end = clock();
    float duration = (float)(end - start);
    cout<<duration<<endl;
    freearray(ppp);
    freearray(dataSetArray);
    freearray(distance);
    freearray1(means);
    
    return 0;
}
