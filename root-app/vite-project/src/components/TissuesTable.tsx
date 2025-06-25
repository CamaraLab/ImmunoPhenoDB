import React, { useMemo, useState, useEffect } from 'react'
import {
    MaterialReactTable,
    type MRT_ColumnDef,
  } from 'material-react-table';
import _ from 'underscore';
import { default as axios } from "axios";
import { select } from 'underscore';

/* 
The TissuesTable.tsx component consists of one Material React Table, with two columns called
Tissue ID  and TissueName.

An API call is made to the /api/tissues endpoint to retrieve all the current tissues in the database.
This JSON is converted into an array object for the MaterialReactTable component.

At least one tissue must be selected in order for downstream API calls to be made.
*/

//  const URL = "http://127.0.0.1:8080"
// const URL = "http://host.docker.internal" // local development, with docker

// Cloud development
// const URL = "http://node:8080"
// const URL = "http://ec2-100-28-86-22.compute-1.amazonaws.com"
// const URL = ""

// Cloud deployment -- new approach (04/17/2025) using react bundled in NGINX
const URL = "" // it should be empty so NGINX intercepts it


// Interfaces
interface TissuesTable {
  idBTO: string;
  tissueName: string;
  id?: string;
};

interface RowSelection {
  [key: string]: boolean;
}

// Helper functions to process tissue data
function convertJSONToArrayObject(dictionary: TissuesTable) {
  const resultspt1: any = []
  for (const [key, values] of Object.entries(dictionary)) {
      for (const [index, value] of Object.entries(values as any[])) {
          (resultspt1[index]  ??= {})[key] = value
      }
  }
  return resultspt1
} 

function includeIDInJSON(transposeCol: any, keyName: string) {
  const resultspt2: any = []
  for (const [key, values] of Object.entries(transposeCol)) {
      let temp: any = values
      temp[keyName] = key
      resultspt2.push(temp)
  }
  return resultspt2
}

const TissuesTable = ({setTissueRows}: {setTissueRows: any}) => {
  const [tissues, setTissues] = useState<TissuesTable[]>([])
  const [rowSelection, setRowSelection] = useState<RowSelection>({});
  const [isError, setIsError] = useState(false)

  // Layout of the MRT Table for Tissues
  const tissuesColumns = useMemo<MRT_ColumnDef<TissuesTable>[]>(
    () => [
        {
            accessorKey: 'idBTO',
            header: 'Tissue ID',
            size: 200
        },
        {
            accessorKey: 'tissueName',
            header: 'Tissue Name',
            size: 300
        }
    ],
    [],
  );

  const fetchTissues =  async () => {
    // Make API call to tissues
    const tissuesRes = await axios.get(`${URL}/api/tissues`)
    
    // Extract JSON result
    const tissuesJSON = tissuesRes.data

    // Convert JSON to Array
    const tissuesArray = convertJSONToArrayObject(tissuesJSON)

    // Add "id" field in Array object
    const tissuesArrayID = includeIDInJSON(tissuesArray, "id")

    // Select all the tissues initially
    let initialSelectedTissues: any = {};
    for (let i=0; i<tissuesArrayID.length; i++){
        let id = tissuesArrayID[i].id
        initialSelectedTissues[id] = true;
    }

    setTissues(tissuesArrayID)              // This is a JSON of all the tissues {0: {idBTO:"", tissueName:""}, 1: {}...}
    setRowSelection(initialSelectedTissues) // This is a JSON of {0: boolean, 1: boolean...}
  }

  useEffect(() => {
      fetchTissues()
  }, []) // empty array means it only runs once
  // source: https://stackoverflow.com/questions/59914158/how-to-do-request-on-page-load-react-hooks

  useEffect (() => {
    // Constantly update selected tissues whenever row state (rowSelection) changes
    const tissueSelectedKeys = Object.keys(rowSelection)
    const selectedTissues = _.map(tissueSelectedKeys, function (item: any) {
        return tissues[item]
      })

    setTissueRows(selectedTissues)

    // If no tissues are selected, raise an error
    if (Object.keys(rowSelection).length === 0) {
        setIsError(true)
    } else {
        setIsError(false)
    }
  }, [rowSelection])

  return (
    <div style={{backgroundColor: 'white'}}>
        <MaterialReactTable
            data={tissues}
            columns={tissuesColumns}
            enableRowSelection
            enableColumnResizing
            onRowSelectionChange={setRowSelection}
            initialState={{pagination: { pageSize: 10, pageIndex: 0}}} // make tissues view 10 results at a time
            muiToolbarAlertBannerProps={
                isError
                  ? {
                      color: 'error',
                      children: 'At least one tissue must be selected',
                    }
                  : undefined
              }
            state={{
                rowSelection,
                showAlertBanner: isError
            }}
        />
    </div>
  )
}

export default TissuesTable