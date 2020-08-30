#include "GLC30_Based.h"


int main()
{
	clock_t start, finish;
	double duration = 0;
	start = clock();
	float startTime = omp_get_wtime();
	GDALAllRegister();

	int GLC30[12] = { 0 };

	int argc = 1;
	//string shp_ID = "n45_35";
	//string tif_ID = "n45_35_2000lc030";
	//string file = "split\\";
	//string shp_head = "E:\\program\\GLC30andDGGRID\\data\\SHP\\";
	//string tif_head = "E:\\program\\GLC30andDGGRID\\data\\TIF\\";
	//string shp_head_ = shp_head + shp_ID + "\\" + file + shp_ID + "_";
	//string tif_head_ = tif_head + tif_ID + "\\" + file + tif_ID + "_";
	char* token = "dggrid_operation";
	char* remainder = "GENERATE_GRID";
	omp_set_num_threads(5);//经比较，线程数定为7是最快的,台式机

//#pragma omp parallel
	{ 
		ofstream outfile;
		int *tmp;
		DgGridPList plist;
		setParam_(&plist, NLEVEL);
//#pragma omp for
		int k = 0;
		ifstream in("D:\\data\\SAMPLE\\GlobeLand30\\SHP\\MRBGRID\\filename.txt");//以TIF中的文件为标准
		string filename;
		string line;
		string tifPath, shpPath, txtPath, outshpPath;
		while (getline(in, line))
		{
			//从txt中读取路径
			tifPath = line;
			cout << line << endl;
			line  = "n45_35_class_10";
			//GetPath(&shpPath, &txtPath,line);
			//string tifPath = "D:\\data\\TIF\\N45_35_2000LC030\\split\\n45_35_2000lc030_" + to_string(k) + ".tif";
			//string shpPath = "D:\\data\\SHP\\N45_35\\split\\n45_35_" + to_string(k) + ".shp";
			//string txtPath = "D:\\data\\TIF\\N45_35_2000LC030\\split\\n45_35_2000lc030_" + to_string(k) + "D.txt";//ISEA4T ISEA4D ISEA4H
			//shpPath = "D:\\data\\SAMPLE\\InputSHP_\\" + line + ".shp";
			//shpPath = "D:\\data\\SAMPLE\\GlobeLand30\\SHP\\intersect\\19736_inter.shp";
			shpPath = "D:\\data\\SAMPLE\\GlobeLand30\\AreaSample\\SHP\\n47-40\\n47_40_9000_latlon.shp";
			//shpPath = "D:\\data\\SAMPLE\\GlobeLand30\\SHP\\input\\" + line + ".shp";
			outshpPath = "D:\\data\\SAMPLE\\GlobeLand30\\Export_Output_edge.shp";
			//txtPath = "D:\\data\\SAMPLE\\GlobeLand30\\SHP\\MRBGRID\\MRBGRIDMAT\\MRBGRID_" +line + ".txt";
			txtPath = "D:\\data\\SAMPLE\\GlobeLand30\\SHP\\MRBGRID\\MRBGRID\\MRBGRID\\MRBGRID_leave_1.txt";//OURSGRID_" + line + ".txt";
																					/*tifPath = "D:\\data\\SAMPLE\\TIF\\test1.tif";*/
			string Snyder_txtPath = "D:\\data\\SAMPLE\\GlobeLand30\\AreaSample\\TXT\\n47-40\\n47_40_9000_SnyderFwd.txt";
			plist.setParam("clip_region_files", (char *)(shpPath.data()));
			//outfile.open(txtPath.data());
			//仅用于计算边界格元
			shpPath = "D:\\data\\tmp\\MBR_n45_35_2000.shp";
			//shpPath = "D:\\data\\SAMPLE\\GlobeLand30\\SHP\\originData\\n45_35_4raster2shp.shp";
			//outshpPath = "D:\\data\\SAMPLE\\GlobeLand30\\Paper2\\n47-35\\EdgeGrid\\n47_35_2010_edgegrid_18.shp";
			outshpPath = "D:\\data\\tmp\\" + to_string(_TYPE_) + "_";
			//outshpPath = "D:\\data\\SAMPLE\\GlobeLand30\\SHP\\n45-35\\EdgeGrid\\n45_35_2010_edgegrid_18.shp";
			string outtxtPath = "D:\\data\\tmp\\" + to_string(_TYPE_) + "_";
			//string outtxtPath = "D:\\data\\SAMPLE\\GlobeLand30\\SHP\\n45-35\\Property\\n45_35_2010_";
			start = clock();
			DGGRID_GLC30(plist, (char*)(shpPath.data()), &tmp, line, (char *)(txtPath.data()), (char *)(outshpPath.data()), (char *)(outtxtPath.data()));
			finish = clock();
			duration = finish - start;
			//仅用于计算snyder投影
			/*for (int i = 24990; i <= 33000; i = i + 30) 
			{
				cout << i << endl;
				shpPath = "D:\\data\\SAMPLE\\GlobeLand30\\AreaSample\\SHP\\n47-35\\n47_35_" + to_string(i) + "_latlon.shp";
				Snyder_txtPath = "D:\\data\\SAMPLE\\GlobeLand30\\AreaSample\\TXT\\n47-35\\n47_35_" + to_string(i) + "_SnyderFwd.txt";
				DGGRID_GLC30(plist, (char*)(shpPath.data()), &tmp, line, (char *)(Snyder_txtPath.data()), (char *)(outshpPath.data()), (char *)(outtxtPath.data()));
			}*/

			
			//////
			//break; 
			//////
			finish = clock();
			duration = finish - start;
//#pragma omp critical
			/*for (int i = 0; i < 12; i++)
			{
				GLC30[i] = GLC30[i] + tmp[i];
				outfile << tmp[i] << endl;
				cout<< tmp[i] << endl;
			}
			cout << k << endl;
			outfile.close();*/
		}
	}
	GDALDestroyDriverManager();
	double endTime = omp_get_wtime();
	double time_tmp;
	time_tmp = endTime - startTime;
	system("pause");
}