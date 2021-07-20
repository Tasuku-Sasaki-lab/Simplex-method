#include <stdio.h>
#include <math.h>
#define epsiron 0.00005
#define ROW 3 // 制約条件の数(基底変数の数) 
#define SLACK 3 ///* SLACK変数の数 */
#define COLUMN 2 ///* 列数(目的関数の変数の数) */
#define M 3 ///* 基底変数の数 */
#define N 5 ///* 全ての変数の数 MとNを使ってプログラムを組んでも構いません */

int main(void)
{
	int i, j, k, end_flag;
	int pivot_column,pivot_row; ///* ピボット変換の基準となる値を入れるための変数 */
	double temp; ///* 値の保持のための変数 */
	int min_val;
	double min_num;

    double A[ROW][COLUMN+SLACK],B[ROW],C[COLUMN+SLACK]; ///* 各係数を入れるための配列の準備 */
    double d[COLUMN+SLACK]; ///* 残差を入れるための配列 */
    int Cbx[ROW]; ///* 基底変数となる番号を入れた配列 */
    double Cb[ROW];// /* 基底変数となる目的関数の係数を入れた配列 */
    FILE *fp; ///* ファイル読み込みのための変数 */
    char fname[30]; ///* ファイル名の入力用変数　*/

	double ob; ///* 目的関数の値の計算用の変数 */

    ///* 係数情報ファイルからの情報取得 */
    printf("Input file name =");
    scanf("%s",fname);

    fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("cannot open file");
        return 0;
    }
    else{
        for(i=0;i<COLUMN;i++){ ///* 目的関数の入力 */
            fscanf(fp,"%lf",&C[i]);
        }
        for(i=0;i<ROW;i++){// /* 各係数の入力 */
            for(j=0;j<COLUMN+SLACK;j++){
                fscanf(fp,"%lf ",&A[i][j]);
            }
            fscanf(fp,"%lf",&B[i]);
        }
        for(i=COLUMN+SLACK-ROW;i<COLUMN+SLACK;i++){ ///* スラック変数部分の目的関数の係数を0にするための処理 */
            C[i]=0.00;
        } 
        for(i=0;i<ROW;i++){ ///* 基底変数の番号の入力 */多分最初は　２，３，４
            fscanf(fp,"%d",&Cbx[i]);
        }
    }

    ///* 初期設定 */
    for(i=0;i<COLUMN+SLACK;i++){ ///* 残差の初期設定 */
        d[i]=C[i];
    }
    for(i=0;i<ROW;i++){ ///* 基底変数の目的関数の係数の入力 */
        Cb[i]=C[Cbx[i]];
    }
	

	end_flag = 0; ///* 全てのdが0以上になったかを判断するためのフラグ */
	while(1){
		i =0;
		
		///* 全てのd（残差）が0以上になったかどうかを調べ，なっていれば結果を出力させるためにend_flagの値を変更する */
		///* 丸め誤差を考慮してepsironを使用する */

		for(j=0;j<COLUMN+SLACK;j++){
			if(d[j] > - epsiron){
				i++;
			}
		}
		if(i == COLUMN+SLACK)
			end_flag=1;



		if (end_flag ==0){ ///* 負の残差がある場合，シンプレックス法を解く */
			min_num = 0;
			min_val = -1;

			////* ピボットする列を決定するため，最も小さいdの値を探索し，ピボット変換する列を決定 */
			for(j=0;j<COLUMN+SLACK;j++){
				if(min_num>d[j]){
					min_num=d[j];
					min_val=j;
				}
			}

			pivot_column = min_val;

			

            ///* ピボット変換するための行を決定 */
			///* nピボットする行は，Aの値が0でなく，B[j]/A[j][pivot_colum]が正の値のうち，最も小さい値 */
			k=0;
			for(j=0;j<ROW;j++){
				if(fabs(A[j][pivot_column])>epsiron && B[j]/A[j][pivot_column]>0){
					if(k==0){
						min_val=j;
						min_num=B[j]/A[j][pivot_column];
						k=1;
					}
					if(B[j]/A[j][pivot_column]<min_num){
						min_num=B[j]/A[j][pivot_column];
						min_val=j;
					}
				}
			}

			pivot_row=min_val;




            //* ガウス・ジョルダン法（掃き出し法）による変換 */
			///* ピボットする行に対する演算 */
			temp=A[pivot_row][pivot_column];
			for(j=0;j<COLUMN+SLACK;j++){
				A[pivot_row][j]/=temp;
			}
			B[pivot_row]/=temp;





			///* それ以外の行に対する演算 */
			///* ピボット変換を行った行以外で演算を行う */
			for(i=0;i<ROW;i++){
				if(i!=pivot_row){
					temp=A[i][pivot_column];
					for(j=0;j<COLUMN+SLACK;j++){
						A[i][j] = A[i][j] - temp*A[pivot_row][j];
					}
					B[i]=B[i]-temp*B[pivot_row];
				}
			}



			///* 残差のところまで忘れずに演算を行うことに注意 */
			temp=d[pivot_column];
			for(j=0;j<COLUMN+SLACK;j++){
				d[j]=d[j]-temp*A[pivot_row][j];
			}



			///* 基底変数の入れ替え */
			Cbx[pivot_row] = pivot_column;
			Cb[pivot_row] = C[Cbx[pivot_row]];


			//printf("%d %d\n",pivot_row,pivot_column);

			for(k=0;k<ROW;k++){ ///* 途中の変換の様子を出力 */
				for(j=0;j<COLUMN+SLACK;j++){
					printf("%lf ",A[k][j]);
				}
				printf("%lf\n",B[k]);
			}

			/*
			for(i=0;i<COLUMN+SLACK;i++){
				printf("%f ",d[i]);
			}
			printf("\n");
			*/


		}
		else{
			///* 全ての残差が0以上になったら結果を出力 */
			ob = 0.000;

			for(k=0;k<ROW;k++){ ///* 目的関数値の計算 */
				ob+= Cb[k] * B[k];
			}
			printf("Objective Function = %d\n",(int)ob);
			for(k=0;k<ROW;k++){ ///* 基底変数の番号とその値の出力 */
				if(Cbx[k] < COLUMN+SLACK-ROW){
					printf("Basic variable: x%d = %lf\n",Cbx[k],B[k]);
				}else{
					printf("Basic variable: Lambda%d = %lf\n",Cbx[k]-(COLUMN+SLACK-ROW),B[k]);
				}
			}

			///* これ以下で総販売金額，総利益など求められている値を出力させる */
			int total=-(int)ob;
			printf("Total Sales = %d yen,",total);
			printf("Profits = %d yen,",(int)(total-200*B[0]-150*B[1]));
			printf("Production volume: P1 = %f gram, P2 = %f gram,",B[0],B[1]);
			printf("Used material volume: M = %f gram,",100*B[1]+180*B[1]);
			printf("Used solvent volume: S = %f litre",1.5*B[0]+B[1]);

			



			break;
		}

	}

	return 0;

}

