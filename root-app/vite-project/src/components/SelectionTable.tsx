import React, { useState, useEffect } from 'react'

import { 
    Button, 
 } from '@mui/material';

import {
    MaterialReactTable,

  } from 'material-react-table';
  
import _ from 'underscore';

/* 
The SelectionTable.tsx component displays the results retrieved from
the RadioSearchField.tsx component into a MaterialReactTable. Here, this is an option to
select rows in the table and then "Add" them using the provided button to
the Cart.tsx component for further use.

This will take in props from RadioSearchField to adjust the column layout in this table.

*/

const SelectionTable = ({searchResults, // This is retrieved from the API call made in RadioSearchField. We are displaying it here
                       searchColumns,
                       searchIsLoading,
                       searchRowSelection,
                       setSearchRowSelection,
                       setCartData,
                       radioValue,
                      
                      //  setSearchColumns,
                      //  setCartColumns
                      }: 
                      {searchResults: any,
                       searchColumns: any,
                       searchIsLoading: boolean,
                       searchRowSelection: any,
                       setSearchRowSelection: any,
                       setCartData: any,
                       radioValue: any,
                      
                      //  setSearchColumns: any,
                      //  setCartColumns: any
                      } ) => {
  // Save which rows are selected
  const [rowSelection, setRowSelection] = useState({});

  // Array of objects we're working with: searchRowSelection
  function addAntibody (antibodyRow: any) {
    const existingAbIDs = searchRowSelection.map((ab: any) => ab.idAntibody)
    // Only add to cart if it doesn't already exist
    if (!existingAbIDs.includes(antibodyRow.idAntibody)) {
        searchRowSelection.push(antibodyRow);
        setSearchRowSelection([...searchRowSelection])
    }
  }

  function addCellType (cellTypeRow: any) {
    const existingCellTypes = searchRowSelection.map((ct: any) => ct.idCL)
    if (!existingCellTypes.includes(cellTypeRow.idCL)) {
        searchRowSelection.push(cellTypeRow);
        setSearchRowSelection([...searchRowSelection])
    }
  }

  const handleAddButton = (e: any) => {

    const selectedKeys = Object.keys(rowSelection)
    const selectedRows = (_.map(selectedKeys, function (item: any) {
        return searchResults[item]
    }))

    if (radioValue === "Antibodies") {
        for (const ab of selectedRows) {
            addAntibody(ab)
        }
    } else if (radioValue === "Cell Types") {
        for (const ct of selectedRows) {
            addCellType(ct)
        }
    }
}

// Constantly update cart to reflect current selection/deselection in the meantime
  useEffect(() => {
    setCartData(searchRowSelection)
  }, [searchRowSelection])

  // Discard selected rows & shopping cart values when you switch radioValues
  useEffect(() => {
    setSearchRowSelection([])
    setCartData([])
  }, [radioValue])

  // Discard selected rows when you re-submit the query
  useEffect(() => {
    setRowSelection({})
  }, [searchResults])

  
  
  return (
    <div>
    <MaterialReactTable
      columns={searchColumns}
      data={searchResults}
      enableRowSelection
      enableColumnResizing
      initialState={{pagination: { pageSize: 5, pageIndex: 0}}}
      onRowSelectionChange={setRowSelection}
      state={{
        rowSelection: rowSelection,
        isLoading: searchIsLoading
      }}
    />

    <Button 
        sx={{ mt: 1, mr: 1 }} 
        type="submit" 
        variant="contained"
        color="primary" 
        disabled={Object.keys(rowSelection).length === 0 ? true : false} 
        onClick={handleAddButton}>
            Add
    </Button>
    </div>
  );
}

export default SelectionTable