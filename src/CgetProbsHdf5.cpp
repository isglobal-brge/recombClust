#include "include/CgetProbsHdf5.h"


std::vector<BLOCKS> getChunkSelection( std::vector<int>val)
{
    
    
    std::vector<int> distances(val.size());
    std::vector<int> result(val.size());
    
    std::adjacent_difference (val.begin(), val.end(), distances.begin());
    /***/
    std::cout << "using default adjacent_difference: ";
    for (int i=0; i<val.size(); i++) std::cout << distances[i] << ' ';
    std::cout << '\n';
    
    // std::cout << "using default result: ";
    // for (int i=0; i<val.size(); i++) std::cout << result[i] << ' ';
    // std::cout << '\n';
    
     /***/
    
    // Some iterator to the begin and end of a -1 block in a source data
    decltype(std::begin(distances)) endOfBlock{};
    
    auto i1 = std::adjacent_find(distances.begin(), distances.end());
    auto i2 = i1;
    
    std::vector<BLOCKS> selmodel;
    
    int block = 0;
    
    while( i1 != distances.end() )
    {
        
        endOfBlock = std::adjacent_find(i1, std::end(distances), std::not_equal_to<int>());
        
        //Push back new blocks taking in to account blocks with single element.
        if( selmodel.empty() && std::distance(distances.begin(), i1) != 0) {

            if( std::distance(distances.begin(), i1) == 1 && distances[0]<=1 ) {
                // add a block for every single element
                for(int i = 0; i<std::distance(distances.begin(), i1); i++)
                {
                    Rcpp::Rcout<<"\n\n DISTANCIA ... : "<<std::distance(distances.begin(), i1)<<"\n\n";
                    selmodel.push_back(BLOCKS());
                    selmodel[block].start = i;
                    selmodel[block].end = i;
                    selmodel[block].len = 1;
                    selmodel[block].startSelection = val.at(selmodel[block].start);
                    selmodel[block].endESelection = val.at(selmodel[block].end);
                    block++;
                }
            } else {
                for(int i = 1; i<std::distance(distances.begin(), i1); i++)
                {
                    selmodel.push_back(BLOCKS());
                    selmodel[block].start = i;
                    selmodel[block].end = i;
                    selmodel[block].len = 1;
                    selmodel[block].startSelection = val.at(selmodel[block].start);
                    selmodel[block].endESelection = val.at(selmodel[block].end);
                    block++;
                };
            }
            
        } else if (!selmodel.empty() && (selmodel[block-1].end + 1 != std::distance(distances.begin(), i1)) ) {
            for(int i = selmodel[block-1].end + 1; i<std::distance(distances.begin(), i1)-1; i++)
            {
                selmodel.push_back(BLOCKS());
                selmodel[block].start = i;
                selmodel[block].end = i;
                selmodel[block].len = 1;
                selmodel[block].startSelection = val.at(selmodel[block].start);
                selmodel[block].endESelection = val.at(selmodel[block].end);
                block++;
            }
        }
        
        
        
        selmodel.push_back(BLOCKS());
        
        // Take in to account where starts sequence
        if(std::distance(distances.begin(), i1) == 0 ) {
            selmodel[block].start = std::distance(distances.begin(), i1);
        } else {
            selmodel[block].start = std::distance(distances.begin(), i1) - 1;
        }
        
        selmodel[block].end = std::distance(distances.begin(), i1) + std::distance(i1, endOfBlock);
        if(selmodel[block].end == val.size()) {
            selmodel[block].end = selmodel[block].end - 1;
        }
        
        selmodel[block].len = selmodel[block].end - selmodel[block].start;
        selmodel[block].startSelection = val.at(selmodel[block].start);
        selmodel[block].endESelection = val.at(selmodel[block].end);
        
        i1 = endOfBlock;
        if(endOfBlock != distances.end()) {
            i1 = std::adjacent_find(endOfBlock++, distances.end());
            i2 = std::find(endOfBlock++, distances.end(), 1);
            if (std::distance(i1, endOfBlock) != std::distance(i2, endOfBlock)) {
                i1 = i2;
            }
        }
        
        block++;
        
        
    }
    
    // Show result to user. Number of blocks and number of elements in each block
    for (unsigned int i{}; i < block; ++i) {
        
        std::cout << "selects[i].start " << selmodel[i].start << " selects[i].end " << selmodel[i].end << " \n";
        std::cout << "Start line : " << selmodel[i].startSelection << "  -  End line :  " << selmodel[i].endESelection << " \n";
        
    }
    
    return(selmodel);
}



//' Compute cluster Recomb freq by mean of voting
//'
//' @param filename with initial probabilities
//' @param group String, inside we have one dataset for each chromosome
//' @param dataset String, chromosome dataset inside the group
//' @param selection overlap ranges
//' @param nCols Integer, number of columns in selection.
//' @param outgroup String output group
//' @param outdataset String output dataset
//' @return NumericVector with recombination probabilites
//' @export
// [[Rcpp::export]]
Rcpp::RObject getProbs_hdf5( std::string filename, 
                             std::string group, 
                             std::string dataset, 
                             Rcpp::RObject selection, 
                             int nCols,
                             Rcpp::Nullable<std::string> outgroup = R_NilValue,
                             Rcpp::Nullable<std::string> outdataset = R_NilValue )
{
    
    std::string strOutGroup,
                strOutDataset;
    IntegerVector dims_out;
    
    H5File* file;
    
    DataSet* pdatasetin = nullptr;
    //..// DataSet* pdatasetout = nullptr;
    DataSet* unlimDataset = nullptr;
    
    try {
        
        H5::Exception::dontPrint();  
        
        // Eigen::Map<Eigen::VectorXd> vselect = Rcpp::as< Eigen::Map<Eigen::VectorXd> >(selection);
        std::vector<int> vselect = Rcpp::as<std::vector<int>>(selection);
        
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
        std::string strinDatasetName = group + "/" + dataset;
        
        std::sort( vselect.begin(), vselect.end());
        
        std::vector<BLOCKS> fullModel = getChunkSelection(vselect);
        
    /*** M'HE QUEDAT AQUÍ I
    
    
    A PARTIR D'AQUÍ S'HA DE REVISAR TOT ... COM FA EL CÀLCUL? PER FILES O PER COLUMNES ???
     TENIR-HO EN COMPTE A L'HORA DE FER LA LECTURA DE LES DADES I EL COMPUT DELS MODELS !!!!'
    
    ***/
        // GESTIONEM ELS SELECTIONS ...
        
        Eigen::VectorXd means (nCols);
        Eigen::MatrixXd dataSelection (vselect.size(), nCols);
        int xoffset = 0;
        
        Rcpp::Rcout<<"\nDataset to search : "<<strinDatasetName<<"\n";
        if( !exists_HDF5_element_ptr(file, strinDatasetName) ) {
            Rcpp::Rcout<<"Dataset doesn't exists, please, create it before proceed";
            file->close();
            return wrap(3);
        }
        
        
            
        Rcpp::Rcout<<" BLOCKS IN THIS SELECTION : "<<fullModel.size()<<"\n";
        
        pdatasetin =  new DataSet(file->openDataSet(strinDatasetName));
        
        int istartMatrixPos = 0;
         
        for (unsigned int i{}; i < fullModel.size(); ++i) {

            // Get number of rows to read
            int nRows = fullModel[i].end - fullModel[i].start ;
            NumericMatrix readeddata(nRows, nCols);
            
            Rcpp::Rcout<<" Rows to be read : "<<nRows<<"\n";
            Rcpp::Rcout<<" Cols to be read : "<<nCols<<"\n";
            
            Rcpp::Rcout<<" Start to read : "<<fullModel[i].startSelection<<"\n";
            Rcpp::Rcout<<" Start to read - 2 : "<<fullModel[i].startSelection - 1<<"\n";
            
        
            // IntegerVector stride = IntegerVector::create(fullModel[i].startSelection, 0 );
            IntegerVector stride = IntegerVector::create(1, 1 );
            IntegerVector block = IntegerVector::create(1, 1 );
            IntegerVector offset = IntegerVector::create(fullModel[i].startSelection - 1 , 0);
            IntegerVector count = IntegerVector::create(nRows, nCols);

            // Read data from hyperslab in RowMajor and addapt to ColMajor
            read_HDF5_matrix_subset( file, pdatasetin, offset, count, stride, block, REAL(readeddata) );
            Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> tmpBlock (REAL(readeddata), nRows, nCols);
            dataSelection.block( istartMatrixPos, 0 , nRows, nCols ) = tmpBlock; // = read chunk
            
            // update write matrix position
            istartMatrixPos = istartMatrixPos + nRows;
        
        }
        
        pdatasetin->close();
        
        /***** M'HE QUEDAT AQUÍ !!!!  
         * FALTARIA REPASSAR-HO TOT, ARA GENERO LA MITJANA PER MODEL
         * FALTARIA GUARDAR LES DADES AL FITXER I REVISAR SI REALMENT
         * HO HE DE CALCULAR PER FILES O COLUMNES
         * REVISAR-HO FENT UNA EXECUCIÓ DE RECOMBCLUST AMB LES MATEIXES DADES
         * PERÒ LA VERSIÓ EN MEMÒRIA QUE ESTÀ FORÇA MES TESTEJADA !!!! '

        // Steps :
        //  1. Existeix dataset?
        //  2.  Si : Afegir una fila ??
        //      No : Crear-lo
        //  3. Tamany dataset ? 
        //  4. Afegir model a la posició que correspongui.
        
        pdatasetin =  new DataSet(file->openDataSet(stroutDatasetName));
        // write_HDF5_matrix_subset()
        
        
        
        ***/
        
        
        IntegerVector ustride = IntegerVector::create(1, 1 );
        IntegerVector ublock = IntegerVector::create(1, 1 );
        IntegerVector uoffset = IntegerVector::create(0, 0);
        
        
        // // STORE THIS DATA !!!
        // // Rcpp::Rcout<<"\nValors mitjanes rowwise : \n"<< dataSelection.rowwise().mean()<<"\n";
        // Rcpp::Rcout<<"\nValors mitjanes colwise : \n"<< dataSelection.colwise().mean()<<"\n";
        
        bool datasetexists = exists_HDF5_element_ptr( file, stroutDatasetName );
        IntegerVector ucount = IntegerVector::create(1, nCols);
        
        // Rcpp::Rcout<< "Existeix ?? : "<<stroutDatasetName<<"\n";
            
        // if dataset doesn't exist --> Create unlimited dataset
        if( datasetexists == 0 ) {
            create_HDF5_unlimited_matrix_dataset_ptr(file, stroutDatasetName, ucount[0], ucount[1], "numeric");
            Rcpp::Rcout<< "L'acabem de CREAR !!!\n \n";
            dims_out[0] = 0; dims_out[1] = 0;
        }
        
        // // extend dataset (if datasets exists previous execution)
        unlimDataset = new DataSet(file->openDataSet(stroutDatasetName));

        if( datasetexists != 0 ) {
            dims_out = get_HDF5_dataset_size(*unlimDataset);
            extend_HDF5_matrix_subset_ptr(file, unlimDataset, 1, 0);
        }
        
        uoffset[0] = dims_out[0];
        
        Rcpp::Rcout<<"Escriurem a : "<<uoffset[0]<< "  -  " <<uoffset[1]<<"\n";
        Rcpp::Rcout<<"Mida a : "<<ucount[0]<< "  -  " <<ucount[1]<<"\n I referent al vector : "<< dataSelection.colwise().mean().size()<<"\n";
        
        Eigen::RowVectorXd toWrite = dataSelection.colwise().mean();
        
        write_HDF5_matrix_subset(file, unlimDataset, uoffset, ucount, ustride, ublock, Rcpp::wrap(toWrite) );  
        
        
        
        
        /** REVISAR AQUÍ QUE ÉS EL QUE HE D'ESCRIURE !!!
         
        FI REVISAR AQUÍ QUE ÉS EL QUE HE D'ESCRIURE !!! **/
        
        
        // Rcpp::Rcout<<"\nAra retornarem... \n";
        

    }catch( FileIException& error ) { // catch failure caused by the H5File operations
        pdatasetin->close();
        unlimDataset->close();
        file->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (File IException)" );
        return wrap(-1);
    } catch( DataSetIException& error ) { // catch failure caused by the DataSet operations
        pdatasetin->close();
        unlimDataset->close();
        file->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataSet IException)" );
        return wrap(-1);
    } catch( DataSpaceIException& error ) { // catch failure caused by the DataSpace operations
        pdatasetin->close();
        unlimDataset->close();
        file->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataSpace IException)" );
        return wrap(-1);
    } catch( DataTypeIException& error ) { // catch failure caused by the DataSpace operations
        pdatasetin->close();
        unlimDataset->close();
        file->close();
        ::Rf_error( "c++ exception Normalize_hdf5 (DataType IException)" );
        return wrap(-1);
    }catch(std::exception &ex) {
        pdatasetin->close();
        unlimDataset->close();
        file->close();
        Rcpp::Rcout<< ex.what();
        return wrap(-1);
    }
    
    
    unlimDataset->close();
    file->close();
    
    
    return(0);
    
    
    
}
