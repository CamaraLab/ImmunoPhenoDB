import React, { useCallback, useMemo, useState, useEffect, useRef } from 'react'
import {
    MaterialReactTable,
} from 'material-react-table';

import { 
    Button, 
    Box,
    Typography,
    Grid,
} from '@mui/material';

import { ExportToCsv } from 'export-to-csv';
import FileDownloadIcon from '@mui/icons-material/FileDownload';

import ReactGA from "react-ga4";

// Modulus function
const mod = (n: number, m: number) => ((n % m) + m) % m;

// For plotly below, we need to unpack data from our object first
function unpack(rows:any, key:any) {
  return rows.map(function(row:any) {return row[key]})
}

/* This needs to handle different types of result tables 
If the radioValue is antibodies, the input will be in the form of a dictionary for both the table and plot
    Accessing these require modelAbKeys, an array which will contain the keys to the dictionary
    We also need the CURRENT key that we want to use to view the table and plot. This is retrieved from tableIndex


*/
const ModelResultsTable = ({radioValue,
                            modelColumns,
                            modelIsLoading,

                            emptyResultsError,
                            setEmptyResultsError, // Set value if table results were empty/bad

                            modelResults,   // For "Cell Types" radioValue
                            
                            
                            modelAbKeys,    // For "Antibodies" radioValue
                            tableIndex,     // For "Antibodies" radioValue
                            setTableIndex,   // For "Antibodies" radioValue
                            modelResultsDict
                           }:
                           {radioValue: any,
                            modelColumns: any,
                            modelResults: any,
                            modelIsLoading: any,
                            emptyResultsError: any,
                            setEmptyResultsError: any,
                            modelAbKeys: any,
                            tableIndex: number,
                            setTableIndex: any,
                            modelResultsDict: any
                           }
                            ) => {

  // CSV properties
  const time = new Date().toISOString();
  const csvOptions = {
    fieldSeparator: ',',
    filename: `${time}`,
    quoteStrings: '"',
    decimalSeparator: '.',
    showLabels: true,
    useBom: true,
    useKeysAsHeaders: true,
    useTextFile: false,
    };
  
  const csvExporter = new ExportToCsv(csvOptions);
  const textFileOptions = {
    fieldSeparator: ',',
    filename: `${time}`,
    quoteStrings: '"',
    decimalSeparator: '.',
    showLabels: true,
    useBom: true,
    useKeysAsHeaders: true,
    useTextFile: true,
    };
  const textFileExporter = new ExportToCsv(textFileOptions)

//   const currentAntibody = modelAbKeys[tableIndex]

  const handleExportDataCSV = () => {
    // Send the event to Google Analytics
    ReactGA.event({
      category: "data_export",
      action: "export_csv_results",
    });

    // If modelAbKeys is empty, then we are dealing with a single 'Cell Types' array
    if (modelAbKeys.length === 0 ) {
        csvExporter.generateCsv(modelResults);
    // If it isn't empty, we're dealing with an 'Antibodies' dictionary
    } else if (modelAbKeys.length >= 1) {
        const currentAntibody = modelAbKeys[tableIndex]
        csvExporter.generateCsv(modelResultsDict[currentAntibody]);
    } else {
        return;
    }
  }

  const handleExportDataTextFile = () => {
    ReactGA.event({
      category: "data_export",
      action: "export_txt_results",
    });

    // If modelAbKeys is empty, then we are dealing with a single 'Cell Types' array
    if (modelAbKeys.length === 0 ) {
      textFileExporter.generateCsv(modelResults);
    // If it isn't empty, we're dealing with an 'Antibodies' dictionary
    } else if (modelAbKeys.length >= 1) {
      const currentAntibody = modelAbKeys[tableIndex]
      textFileExporter.generateCsv(modelResultsDict[currentAntibody]);
    } else {
        return;
    }
  }

  // Event handlers for forward and backward buttons.
  // Used when displaying plots for "Antibodies" radioValue
  const forwards = useCallback(() => 
    setTableIndex((state: number) => mod(state - 1, modelAbKeys.length))
  , [setTableIndex, modelAbKeys])
  
  const backwards = useCallback(() => 
    setTableIndex((state: number) => mod(state + 1, modelAbKeys.length))
  , [setTableIndex, modelAbKeys])

  const handleForwards = (e: any) => {
    e.preventDefault()
    if (radioValue === "Antibodies") {
        forwards()
    } 
  }

  const handleBackwards = (e: any) => {
      e.preventDefault()
      if (radioValue === "Antibodies") {
          backwards()
      } 
  }


  // Prepare data to display in the table, depending on the radioValue chosen
  let modelData: any = []
  if (radioValue == "Antibodies") {
    if (modelAbKeys.length > 0) {
        const currentAntibody = modelAbKeys[tableIndex]
        console.log("We're using antibody data. Current antibody:", currentAntibody)
        modelData = modelResultsDict[currentAntibody]
        if (modelData == null) {
            modelData = [] // If there were no results to begin with, but we kept modelAbKeys to display an error - set results to empty array
        }
    }
  } else if (radioValue == "Cell Types") {
        console.log("We're using cell type data")
        modelData = modelResults
  }
  console.log("modelData for table:", modelData)


  // Check if table is empty only AFTER Cart.tsx finishes getting data for table
  useEffect(() => {
    // For "Antibodies" radioValue. Don't show error on intial load, so check if modelAbKeys has more than one antibody key
    // At the same time, check that no results exist for "Cell Types"
    // Finally, check that modelData extracted is also empty
    if (modelIsLoading == false && (modelAbKeys.length >= 1 && radioValue === "Antibodies") && modelData.length === 0) {
      console.log("trigger empty results error for antibodies")
      setEmptyResultsError(true)
    } 
    // For "cell types". Check that no results exist for "Antibodies" while checking modelResults to be empty too
    // Finally check that modelData extracted is also empty
    else if (modelIsLoading == false && (radioValue === "Cell Types" && modelResults.length == 0) && modelData.length === 0){
      console.log("trigger empty results error for cell types")
      setEmptyResultsError(true)
    }
    else {
      console.log("No empty results to trigger here!")
      setEmptyResultsError(false)
    }
  }, [modelAbKeys, modelData])


  useEffect(() => {

  })



  console.log("modelIsLoading:", modelIsLoading)
  console.log("modelAbKeys.length:", modelAbKeys.length)
  console.log("modelResults.length:", modelResults.length)
  console.log("modelData.length:", modelData.length)
  console.log("emptyResultsError:", emptyResultsError)

  useEffect(() => {
    setEmptyResultsError(false)
  }, [radioValue])

  let error_message: any = {}
  if (radioValue == "Antibodies" && emptyResultsError == true) {
    error_message =  
        {color: 'error',
         children: `No results found for antibody ${modelAbKeys[tableIndex]}`}
  } else if (radioValue == "Cell Types" && emptyResultsError == true) {
    error_message = {
      color: 'error', 
      children: 'No results found'
    }
  } 
  
  return (
    <div>
        <div>
            <MaterialReactTable
                columns={modelColumns}
                data={modelData}
                state={{
                    isLoading: modelIsLoading,
                    showAlertBanner: emptyResultsError
                }}
                initialState={{pagination: { pageSize: 5, pageIndex: 0}}}
                enableColumnResizing
                muiToolbarAlertBannerProps={
                  (emptyResultsError)
                    ? error_message
                    : {}
                    }

                renderTopToolbarCustomActions={({ table }) => (
                  <div>
                    <Box 
                      display="flex" 
                      justifyContent="space-between" 
                      alignItems="center" 
                      padding={1}
                      gap={2}
                    >
                      <Box>
                        <Button
                          color="primary"
                          //export all data that is currently in the table (ignore pagination, sorting, filtering, etc.)
                          onClick={handleExportDataCSV}
                          startIcon={<FileDownloadIcon />}
                          variant="contained"
                          disabled={modelIsLoading || emptyResultsError || !modelData.length}
                        >
                          Export As CSV
                        </Button>
              
                        <Button
                          color="secondary"
                          //export all data that is currently in the table (ignore pagination, sorting, filtering, etc.)
                          onClick={handleExportDataTextFile}
                          startIcon={<FileDownloadIcon />}
                          variant="contained"
                          disabled={modelIsLoading || emptyResultsError || !modelData.length}
                          style={{ marginLeft: '10px' }}
                        >
                          Export As TXT
                        </Button>
                      </Box>
              
                      <Box flexGrow={1} textAlign="center">
                        <b>
                          {!emptyResultsError && modelAbKeys.length >= 1 && radioValue === "Antibodies"
                            ? "Viewing: " + modelAbKeys[tableIndex]
                            : ""}
                        </b>
                      </Box>
                    </Box>
                  </div>
                  )}
                />

            <br></br>

            {radioValue === "Antibodies" ?
            <div>
              <Button  
                  sx={{ mt: 1, mr: 1 }} 
                  onClick={handleForwards}
                  variant="contained"
                  disabled={modelAbKeys.length <= 1}>
                      Backward
              </Button>
              <Button 
                  sx={{ mt: 1, mr: 1 }} 
                  onClick={handleBackwards}
                  variant="contained"
                  disabled={modelAbKeys.length <= 1}>
                      Forward
              </Button>
            </div> : <b></b>}
        </div>
    </div>
  )
}

export default ModelResultsTable