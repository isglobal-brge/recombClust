#include "include/CWriteJson.h"


using namespace Rcpp;
using namespace std;


// Write json configuration file for HiCBricks HDF5
int write_jsonBrick_file(std::string strfilename, RObject params)
{
   try
   {
      
      int iresolution = as<List>(params)["resolution"];
      int ilengths = as<List>(params)["length"];
      int idimensions = as<List>(params)["dimensions"];
      std::string strpath = as<List>(params)["path"];
      std::string strchromosomes = as<List>(params)["chromosomes"];
      std::string strhdf5file = as<List>(params)["HDF5File"];
      

      std::ofstream jsonfile (strfilename);
      
      if (jsonfile.is_open())
      {
         jsonfile << "{\n";
         jsonfile << "\t\"headers\": {\n";
         jsonfile << "\t\t\"file_prefix\": [\n";
         jsonfile << "\t\t\t\"invSim\"\n";
         jsonfile << "\t\t],\n";
         
         jsonfile << "\t\t\"project_directory\": [\n";
         jsonfile << "\t\t\t\"" + strpath.substr(0, strpath.size()-1) + "\"\n";
         jsonfile << "\t\t],\n";
         // M'he quedat on s'ha de posar el directori....
         jsonfile << "\t\t\"experiment_name\": [\n";
         jsonfile << "\t\t\t\"Inversion Simulation\"\n";
         jsonfile << "\t\t],\n";
         
         jsonfile << "\t\t\"resolutions\": [\n";
         jsonfile << "\t\t\t\"" +   std::to_string(iresolution)  + "\"\n";
         jsonfile << "\t\t],\n";
         
         jsonfile << "\t\t\"chromosomes\": [\n";
         jsonfile << "\t\t\t\"" + strchromosomes + "\"\n";
         jsonfile << "\t\t],\n";
         
         jsonfile << "\t\t\"lengths\": [\n";
         jsonfile << "\t\t\t" + std::to_string(ilengths)  + "\n";
         jsonfile << "\t\t],\n";
         
         jsonfile << "\t\t\"queues\": [\n";
         jsonfile << "\t\t\tnull\n";
         jsonfile << "\t\t]\n";
         
         jsonfile << "\t},\n";
         
         
         jsonfile << "\t\"matrix_info\": {\n";
         jsonfile << "\t\t\"" + strchromosomes + "_" + strchromosomes + "_" + std::to_string(iresolution) + "\": {\n";
         
         jsonfile << "\t\t\t\"chrom1\":[\n";
         jsonfile << "\t\t\t\t\""+strchromosomes+"\"\n";
         jsonfile << "\t\t\t],\n";
         
         jsonfile << "\t\t\t\"chrom2\":[\n";
         jsonfile << "\t\t\t\t\""+strchromosomes+"\"\n";
         jsonfile << "\t\t\t],\n";
         
         jsonfile << "\t\t\t\"resolution\":[\n";
         jsonfile << "\t\t\t\t\"" + std::to_string(iresolution) + "\"\n";
         jsonfile << "\t\t\t],\n";
         
         jsonfile << "\t\t\t\"dimensions\":[\n";
         jsonfile << "\t\t\t\t" + std::to_string(idimensions) + ",\n";
         jsonfile << "\t\t\t\t" + std::to_string(idimensions) + "\n";
         jsonfile << "\t\t\t],\n";
         
         jsonfile << "\t\t\t\"lengths\":[\n";
         jsonfile << "\t\t\t\t" + std::to_string(ilengths) + ",\n";
         jsonfile << "\t\t\t\t" + std::to_string(ilengths) + "\n";
         jsonfile << "\t\t\t],\n";
         
         jsonfile << "\t\t\t\"mat_type\":[\n";
         jsonfile << "\t\t\t\t\"cis\"\n";
         jsonfile << "\t\t\t],\n";
         
         jsonfile << "\t\t\t\"filename\":[\n";
         // jsonfile << "\t\t\t\t\"invSim_" + std::to_string(iresolution) + "_"+strchromosomes+"_vs_"+strchromosomes+".brick\"\n";
         jsonfile << "\t\t\t\t\"" + strhdf5file + "\"\n";
         jsonfile << "\t\t\t]\n";
         
         jsonfile << "\t\t}\n";
         jsonfile << "\t}\n";
         jsonfile << "}\n";
         jsonfile.close();
      }
      
      else Rcpp::Rcout << "Unable to open file";
      
   }
   catch(std::exception const& e){
      Rcpp::Rcout << "There was an error: " << e.what() << std::endl;
   }
   return 0;
   
}



/*** R

*/
