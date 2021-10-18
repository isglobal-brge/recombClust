#include "include/CgetProbsHdf5.h"

using namespace Rcpp;

//' Compute cluster Recomb freq by mean of voting
//'
//' @param mat with initial probabilities
//' @param sel overlap ranges
//' @return NumericVector with recombination probabilites
//' @export
// [[Rcpp::export]]
Rcpp::RObject getProbs_hdf5( std::string filename, 
                             std::string group, 
                             Rcpp::RObject selection, 
                             Rcpp::Nullable<std::string> outgroup = R_NilValue,
                             Rcpp::Nullable<std::string> outdataset = R_NilValue )
{
    
    std::string strOutGroup,
                strOutDataset;
    
    H5File* file;
    
    DataSet* pdatasetin;
    DataSet* pdatasetout;
    
    try {
        
        H5::Exception::dontPrint();  
        
        if(!ResFileExist_filestream(filename)) {
            Rcpp::Rcout<<"\nERROR - File doesn't exits\n";  
            return Rcpp::wrap(3);
        }
        
        file = new H5File( filename, H5F_ACC_RDWR );
        
        if(outgroup.isNull()) { strOutGroup = group; } 
        else {   strOutGroup = Rcpp::as<std::string>(outgroup); }
        
        if(outdataset.isNull()) { strOutDataset = "FinalModels"; } 
        else {   strOutDataset = Rcpp::as<std::string>(outdataset); }
        
        // Merge all possible chromosomes results in one file to get complete
        // data?? or get one Result dataset for each chromosome data???
        
        std::string stroutDatasetName = strOutGroup + "/" + strOutDataset;
        std::string strinDatasetName = group ;
        
        
        
        
       
        
        // Eigen::MatrixXd dmat = Rcpp::as<Eigen::MatrixXd>(mat);
        // Eigen::VectorXd dsel = Rcpp::as<Eigen::VectorXd>(sel);
        // 
        // Eigen::MatrixXd dmatsel = Eigen::MatrixXd::Zero(dmat.rows(), dsel.size());
        // 
        // for( int i=0; i<dsel.size(); i++){
        //     dmatsel.col(i) = dmat.col(dsel[i]-1);
        // }
        // 
        // return(wrap(dmatsel.rowwise().mean()));
        
        
        
        
        
    }catch( FileIException& error ) { // catch failure caused by the H5File operations
        file->close();
        pdatasetin->close();
        pdatasetout->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (File IException)" );
        return wrap(-1);
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        file->close();
        pdatasetin->close();
        pdatasetout->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataSet IException)" );
        return wrap(-1);
    } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        file->close();
        pdatasetin->close();
        pdatasetout->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataSpace IException)" );
        return wrap(-1);
    } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        file->close();
        pdatasetin->close();
        pdatasetout->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataType IException)" );
        return wrap(-1);
    }catch(std::exception &ex) {
        file->close();
        pdatasetin->close();
        pdatasetout->close();
        Rcpp::Rcout<< ex.what();
        return wrap(-1);
    }
    
    
    
    
    
    
}
