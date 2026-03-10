import React, { useCallback, useMemo, useState, useEffect } from 'react'

import { 
    Button, 
    Box,
    Tooltip,
    IconButton,
 } from '@mui/material';

import {
    MaterialReactTable,
    type MRT_ColumnDef,
    type MRT_Row,
  } from 'material-react-table';
  
import _ from 'underscore';
import { default as axios } from "axios";
import { Delete } from '@mui/icons-material';
import ReactGA from "react-ga4";

/* 
The Cart.tsx component consists of one Material React Table, with three columns called:
Actions, Antibody ID/Cell Ontology ID, and Antibody Target/Cell Type.

There are two action buttons below the table, called "Clear" and "Run"
"Clear" removes all rows in the table
"Run" calls one of two API calls, depending on values in the table
If antibodies are present, /api/findcelltypesweb and /api/plotcelltypesweb are called
If cell types are present, /api/findabsweb and /api/plotabsweb are called

In the Actions column, there is the action to "Delete"
In the Antibody ID/Cell Ontology ID column, there are the RRIDs of each antibody OR the cell ontology
In the Antibody Target/Cell Type column, there are the antibody targets (eg. CD4) 
 
At least one antibody/cell type must be in the cart for either "Clear" or "Run" to be enabled
*/

// const URL = "http://127.0.0.1:8080"
// const URL = "http://host.docker.internal" // local development, with docker

// Cloud deployment
// const URL = "http://node:8080" doesn't work
// const URL = "http://ec2-100-28-86-22.compute-1.amazonaws.com" doesn't work
// const URL = ""

// Cloud deployment -- new approach (04/17/2025) using react bundled in NGINX
const URL = "" // it should be empty so NGINX intercepts it

// Helper functions 
function convertToArrayObject(dictionary: any) {
  const resultspt1: any = []
  for (const [key, values] of Object.entries(dictionary)) {
      for (const [index, value] of Object.entries(values as any[])) {
          (resultspt1[index]  ??= {})[key] = value
      }
  }
  return resultspt1
} 

function includeID(transposeCol: any, keyName:string) {
  const resultspt2: any = []
  for (const [key, values] of Object.entries(transposeCol)) {
      let temp: any = values
      temp[keyName] = key
      resultspt2.push(temp)
  }
  return resultspt2
}

const Cart = ({cartData,
              cartColumns,
              tissueRows,
              setCartData,
              setSearchRowSelection,
              radioValue,
              setEmptyResultsError,
              setModelColumns,
              setModelResults, // Cell Types
              setModelIsLoading,
              setModelPlotIsLoading,
              setModelPlotResults, // Cell Types
              setModelResultsDict,
              setModelPlotResultsDict,
              setModelAbKeys, 
              }:
              {cartData: any,
              cartColumns: any,
              tissueRows: any,
              setCartData: any,
              setSearchRowSelection: any,
              radioValue: any,
              setEmptyResultsError: any,
              setModelColumns: any,
              setModelResults: any,
              setModelIsLoading: any,
              setModelPlotIsLoading: any,
              setModelPlotResults: any,
              setModelResultsDict: any,
              setModelPlotResultsDict: any,
              setModelAbKeys: any
              }) => {

  // State to hide/show 'Run' button underneath the table
  const [runDisabled, setRunDisabled] = useState(false);

  const AntibodyLMMColumns = [
      {
        accessorKey: 'id',
        header: 'Antibody ID',
        // size: 125,
      },
      {
        accessorKey: 'target',
        header: 'Antibody Target',
      },
      {
        accessorKey: 'coeff',
        header: 'Coefficient',
      },
      {
        accessorKey: 'stderr',
        header: 'Standard Error',
      },
      {
        accessorKey: 'p_val',
        header: 'P value',
      },
      {
        accessorKey: 'q_val',
        header: 'Q value',
      },
      // {
      //   accessorKey: 'other',
      //   header: 'Background Cells',
      // },
      //end
    ]
  
  const CellTypeLMMColumns = [
      {
        accessorKey: 'id',
        header: 'Cell Ontology ID',
        // size: 125,
      },
      {
        accessorKey: 'cellType',
        header: 'Cell Type',
      },
      {
        accessorKey: 'coeff',
        header: 'Coefficient',
      },
      {
        accessorKey: 'stderr',
        header: 'Standard Error',
      },
      {
        accessorKey: 'p_val',
        header: 'P value',
      },
      {
        accessorKey: 'q_val',
        header: 'Q value',
      },
      {
        accessorKey: 'expressed',
        header: 'Expressed',
      },
      //end
    ]

  const format_findabsweb_columns = (AntibodyLMMColumns: any, findabsweb_data: any) => {
    // Make copy of original AntibodyLMMColumn names
    const findabsweb_columns_copy = [...AntibodyLMMColumns]
          
    // Find all keys containing "CL"
    const findabsweb_data_idCL_keys = Object.keys(findabsweb_data)
    const idCL_keys = findabsweb_data_idCL_keys.filter((element: string) => element.includes("CL:"))
    
    // Add in their column names to the end
    for (const key of idCL_keys) {
      const tempColumn = {
        accessorKey: `${key}`,
        header: `${key}`
      }

      const currentLength = findabsweb_columns_copy.length
      
      // Checking if antibody Key already exists
      const isPresent = findabsweb_columns_copy.some((col) => (col.accessorKey === key))
      
      // If it doesn't already exist, add it to the position right before the last element
      if (!isPresent) {
        findabsweb_columns_copy.splice((currentLength - 1), 0, tempColumn)
      } else {
        // If antibody key already exists, don't add and skip
        continue
      }
    }
    return findabsweb_columns_copy
  }

  const reorder_plotabsweb_using_findabsweb = (plotabsweb_data_array: any, findabsweb_data_array_id: any) => {
    // Get the antibodies from the plotabsweb_data and findabsweb. We use the order in findabsweb for plotting
    const plotabsweb_data_antibodies = plotabsweb_data_array.map((item: any) => item.idAntibody)
    const findabsweb_data_antibodies = findabsweb_data_array_id.map((item: any) => `${item.id} (${item.target})`);
    // Edge case: if there's fewer antibodies in the plotabsweb compared to findabsweb, filter them
    const boxplot_x_axis_order = findabsweb_data_antibodies.filter((item: any) => plotabsweb_data_antibodies.includes(item))
    // Re-order the plotabsweb_data_array 
    const plotabsweb_data_array_reordered = plotabsweb_data_array.sort((a: any, b: any) => {
      return boxplot_x_axis_order.indexOf(a.idAntibody) - boxplot_x_axis_order.indexOf(b.idAntibody);
    });
    return plotabsweb_data_array_reordered
  }

  // Perform two API calls:
  // 1. First one runs the LMM and gets the results back into a table
  // 2. Second one gets the plotting statistics needed for the boxplot
  const handleRunLMM = async (e: any) => {
    // Send the event to Google Analytics
    ReactGA.event({
      category: "query",
      action: "perform_linear_mixed_model",
    });

    // Prevent reload on page when clicked
    e.preventDefault()

    // Ensure the "Run" button is disabled during the duration of the API calls below
    setRunDisabled(true)

    // If dealing with cell types in the cart, call /api/findabsweb and /api/plotabsweb endpoints
    if (radioValue === "Cell Types") {
      // console.log("cell types in cart. running findabsweb and plotabsweb")

      // Prepare payload for first API call for ModelResultsTable.tsx
      const findabsweb_idCLs = cartData.map((cartRow: any) => cartRow.idCL)
      const findabsweb_tissues = tissueRows.map((tissueRow: any) => tissueRow.idBTO)
      const findabsweb_payload = {
        "idCL": findabsweb_idCLs,
        "idBTO": findabsweb_tissues
      }

      // Clear existing column configuration and results in ModelResultsTable.tsx
      setModelColumns([])
      setModelResults([])
      setModelPlotResults([])

      setModelResultsDict([])
      setModelPlotResultsDict([])
      setModelAbKeys([])

      // setEmptyResultsError(false)

      try {
        // Start loading animation for ModelResultsTable
        setModelIsLoading(true)
        
        // Perform API call
        const findabsweb_res = await axios.post(`${URL}/api/findabsweb`, findabsweb_payload)
        const findabsweb_data = findabsweb_res.data

        // Successful call will produce a non-string response. A string indicates an error
        if (typeof findabsweb_data !== "string") {
          // 1. Format the columns
          // NOTE: This API call creates a custom-length of column names for ModelResultsTable. Format it here
          
          // Make copy of original AntibodyLMMColumn names
          const findabsweb_columns_copy = format_findabsweb_columns(AntibodyLMMColumns, findabsweb_data)
          
          setModelColumns(findabsweb_columns_copy)
          // console.log("findabsweb_columns_copy:", findabsweb_columns_copy)

          // 2. Format the data
          const findabsweb_data_array = convertToArrayObject(findabsweb_data)
          const findabsweb_data_array_id = includeID(findabsweb_data_array, "id")

          setModelResults(findabsweb_data_array_id)
          // console.log("findabsweb_data_array_id:", findabsweb_data_array_id)

          setModelIsLoading(false)

          // Prepare payload for second API call for ModelResultsPlot.tsx
          const plotabsweb_payload = {
            "abs": findabsweb_data_array_id.map((row: any) => row.id), // Extract all antibody IDs from result
            "idcls": findabsweb_idCLs
          }
          
          // Start loading animation for ModelResultsTable
          setModelPlotIsLoading(true)

          // Perform API call and format data
          const plotabsweb_res = await axios.post(`${URL}/api/plotabsweb`, plotabsweb_payload)
          const plotabsweb_data = plotabsweb_res.data
          const plotabsweb_data_array = convertToArrayObject(plotabsweb_data)
          // setModelPlotResults(plotabsweb_data_array) -- don't set it so soon, we can re-order them

          const plotabsweb_data_array_reordered = reorder_plotabsweb_using_findabsweb(plotabsweb_data_array, findabsweb_data_array_id)
          setModelPlotResults(plotabsweb_data_array_reordered)

          setModelPlotIsLoading(false)
        
        } else if (typeof findabsweb_data === "string") {
          // Error string
          // console.log("No results found!")
          setModelIsLoading(false)
          setModelResults([])
          setModelPlotResults([])
          setModelPlotIsLoading(false)
          // setEmptyResultsError(true)
        }

      } catch (error) {
        setModelIsLoading(false)
        setModelResults([])
        setModelPlotResults([])
        setModelPlotIsLoading(false)
        // setEmptyResultsError(true)
        // console.log("findabsweb/plotabsweb error:", error)
      } finally {
        // Re-enable Run button
        setRunDisabled(false)
      }
    } 
    
    // If dealing with antibodies in cart, call findcelltypesweb and plotcelltypes
    else if (radioValue == "Antibodies") {
      // console.log("antibodies in cart. running findcelltypesweb and plotcelltypesweb")
     
      // Prepare payload for first API call for ModelResultsTable.tsx
      const findcelltypesweb_antibodies = cartData.map((cartRow: any) => cartRow.idAntibody)
      const findcelltypesweb_tissues = tissueRows.map((tissueRow: any) => tissueRow.idBTO)
      const findcelltypesweb_payload = {
        "ab": findcelltypesweb_antibodies,
        "idBTO": findcelltypesweb_tissues
      }

      // Clear existing column configuration and results in ModelResultsTable.tsx
      setModelColumns([])
      setModelResults([])
      setModelPlotResults([])

      setModelResultsDict([])
      setModelPlotResultsDict([])
      setModelAbKeys([])

      // Start loading animation for ModelResultsTable
      setModelIsLoading(true)
      // setEmptyResultsError(false)

      // We will need a way to keep track of the antibodies when it comes to plotting the boxplots
      // Since that will be needed in our payload to plotcelltypesweb
      // setModelAbKeys([])
  
      // Perform API call
      // Ensure the "Run" button is disabled during the duration of the API calls below
      setRunDisabled(true)

      try {
        const findcelltypes_res = await axios.post(`${URL}/api/findcelltypesweb`, findcelltypesweb_payload)
        const findcelltypes_data = findcelltypes_res.data // This response will be a dictionary of string-JSONs

        // Since we do NOT have a single table that holds the results for each antibody (we ran them separately)
        // We will store them in an object (dictionary), so that we remember which antibody has which plotting result
        const antibody_table_dict: any = {}

        // Let's also keep track of which antibodies we received from the server
        // These will be in the form of common antibody names ("AB_123", "AB_234", ...)
        // This will help us later when plotting AND accessing the correct table values
        const antibody_id_key_array = []
        
        // UNLIKE find antibodies, we do NOT return a string message if there are no results
        // It will just be a key['AB_28000923']: value[Empty dictionary]
        // So there will always be a response in a dictionary
        if (Object.keys(findcelltypes_data).length > 0) {
          // Iterate over each key value
          for (const [key, value] of Object.entries(findcelltypes_data)) {
            antibody_id_key_array.push(key) // Store all antibody keys, regardless of their result

            const value_to_parse: any = (value)
            // If value is an empty dictionary, it is NOT a string-JSON, so we CANNOT run JSON.parse() over it
            // So if we get a value that does have a length of 0 (means empty dictionary, not a string-JSON)
            // Then we KNOW that it's an empty dictionary -> no results for this antibody
            if (Object.keys(value_to_parse).length === 0 && typeof(value_to_parse) !== "string") {
              // If we got an empty dictionary, then there were no results for that antibody
              // Add an empty array for that antibody key, so map() does not crash
              antibody_table_dict[key] = []     
            } else if (typeof(value_to_parse) === "string") {
              // Otherwise, we got a successful result, meaning we have a string-JSON
              // Parse, transform, and format it and push it onto our result dictionary of objects
              const tempValue = JSON.parse(value_to_parse)
              const tempValue2 = convertToArrayObject(tempValue)
              const tempValue3 = includeID(tempValue2, "id")

              antibody_table_dict[key] = tempValue3
            }
          }

          // console.log("antibody_table_dict:", antibody_table_dict)
          // console.log("antibody_id_key_array:", antibody_id_key_array)
        }
        // After we are done with our API response, set the columns for cell types
        setModelColumns(CellTypeLMMColumns)

        // Store results for antibodies
        setModelResultsDict(antibody_table_dict)
        setModelAbKeys(antibody_id_key_array)

        // We can set the table loading animation to false
        setModelIsLoading(false)
        

        // Before getting the plotting data for each antibody, check that our result dict is not empty
        // const antibody_plot_dict: any = {}
        const antibody_plot_dict = new Map(); // using Map here

        // We only want to get the plotting result for antibodies that had a valid result (not empty dictionary or array)
        if (Object.keys(antibody_table_dict).length > 0) {
          setModelPlotIsLoading(true)

          // Make individual API calls for each antibody, and store results in antibody_plot_dict
          for (const ab_key of Object.keys(antibody_table_dict)) {
            // For only antibodies that have valid results. If an antibody has an empty array, do NOT make a plot API call
            if (antibody_table_dict[ab_key].length == 0) {
              // antibody_plot_dict[ab_key] = [] // This will be empty
              antibody_plot_dict.set(ab_key, []); // using Map here
              continue // Skip to next antibody
            }

            // Extract all cell type ids for that antibody
            const celltype_id_array = antibody_table_dict[ab_key].map((item: any) => item.id)
            // Construct a payload for plotcelltypesweb
            const plotcelltypesweb_payload = {
              "ab": ab_key,
              "idcls": celltype_id_array
            }
            
            const plotcelltypesweb_res = await axios.post(`${URL}/api/plotcelltypesweb`, plotcelltypesweb_payload)
            const plotcelltypesweb_data = plotcelltypesweb_res.data

            // Check response of plotting API. If it failed, it would return an empty dictionary
            if (Object.keys(plotcelltypesweb_data).length > 0) {
              // Format the data
              const plotcelltypesweb_data_array = convertToArrayObject(plotcelltypesweb_data)
              // console.log("plotcelltypesweb_data_array:", plotcelltypesweb_data_array)

              // Re-order cell types for an antibody based on order in antibody_table_dict for plotting
              const plotcelltypesweb_data_celltypes = plotcelltypesweb_data_array.map((item: any) => item.idCL)
              const findcelltypesweb_data_celltypes = antibody_table_dict[ab_key].map((item: any) => `${item.cellType} (${item.id})`);
              // Edge case: if there's fewer cell types in the plotcelltypesweb compared to findcelltypesweb, filter them
              const boxplot_x_axis_order = findcelltypesweb_data_celltypes.filter((item: any) => plotcelltypesweb_data_celltypes.includes(item))
              // Re-order the plotcelltypesweb_data_array 
              const plotcelltypesweb_data_array_reordered = plotcelltypesweb_data_array.sort((a: any, b: any) => {
                return boxplot_x_axis_order.indexOf(a.idCL) - boxplot_x_axis_order.indexOf(b.idCL);
              });

              // Add re-ordered results to antibody_plot_dict
              // antibody_plot_dict[ab_key] = plotcelltypesweb_data_array_reordered
              antibody_plot_dict.set(ab_key, plotcelltypesweb_data_array_reordered); // using Map here
            }
          }
          setModelPlotResultsDict(antibody_plot_dict)
          setModelPlotIsLoading(false)
          // console.log("antibody_plot_dict:", antibody_plot_dict)

        } else {
          // The result dict was empty, meaning no cell types were found for any antibodies provided
          // For this, we don't want to do anything. Keep the antibody_object_dict as-is with antibody keys and empty arrays
          // We will check later during the plot to see if that antibody contains an empty array
          // If so, then an error will be raised for the Table and the plot will be hidden
          // setEmptyResultsError(true)

          // SEE above: if (antibody_table_dict[ab_key].length == 0)
        }


      } catch (error) {
        // console.log("findcelltypesweb/plotcelltypesweb error:", error)
        setModelIsLoading(false)
        setModelResultsDict([])
        setModelPlotResultsDict([])
        setModelAbKeys(findcelltypesweb_antibodies)
        setModelPlotIsLoading(false)
        // setEmptyResultsError(true)
      } finally {
        // Re-enable run button
        setRunDisabled(false)
      }
      

    }

  }
  
  const handleDeleteRow = useCallback(
    (row: MRT_Row<any>) => {
        //send api delete request here, then refetch or update local table data for re-render
        cartData.splice(row.index, 1);
        setCartData([...cartData]);
    },
    [cartData],
  );

  const handleClearAll = (e: any) => {
    setCartData([])
    setSearchRowSelection([])
  }

  useEffect(() => {
    if ((Object.keys(tissueRows).length === 0) || 
        (Object.keys(tissueRows).length !== 0 && Object.keys(cartData).length === 0)) {
        // Disable the 'Run' button if there are no tissues, and/or there are no items in the cart
        setRunDisabled(true)
    } else {
        setRunDisabled(false)
    }
  }, [tissueRows, cartData])

  return (
    <div>
      <div style={{backgroundColor: 'white'}}>
        <MaterialReactTable
          displayColumnDefOptions={{
            'mrt-row-actions': {
              muiTableHeadCellProps: {
                align: 'center',
              },
              size: 40,
            },
          }}
          columns={cartColumns}
          data={cartData}
          initialState={{pagination: { pageSize: 5, pageIndex: 0}}} // making shopping cart start at 5 results per page
          enableRowActions={true}
          renderRowActions={({ row, table }) => (
            <Box 
                display="flex"
                justifyContent="center"
                alignItems="center">
              <Tooltip arrow placement="right" title="Delete">
                <IconButton color="error" onClick={() => handleDeleteRow(row)}>
                  <Delete />
                </IconButton>
              </Tooltip>
            </Box>
          )}
          // layoutMode='grid'
          />
      </div>
    
    <div>
    <Button 
    sx={{ mt: 1, mr: 1 }} 
    type="submit" 
    variant="contained"
    color="error" 
    onClick={handleClearAll}
    disabled={Object.keys(cartData).length === 0 ? true : false}>
        Clear
    </Button>

    <Button
    sx={{ mt: 1, mr: 1 }} 
    type="submit" 
    variant="contained" 
    color="success"
    onClick={handleRunLMM}
    disabled={runDisabled}>
        Run
    </Button>
    </div>
  </div>
  )
}

export default Cart