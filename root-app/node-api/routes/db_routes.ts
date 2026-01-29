import express, { Express, Router, Request, Response } from 'express'
import {idCLs, 
        tissues,
        experiments,
        findAbs, 
        plotAbs, 
        findCellTypes, 
        plotCellTypes, 
        findExperiments, 
        whichAntibodies,
        whichCellTypes,
        whichExperiments,
        stveaReference,
        convertCellType,
        databaseStatistics,
        antibodiesTable,
        antibodyPanelReference,
        findAbsWeb,
        plotAbsWeb,
        findCellTypesWeb,
        plotCellTypesWeb,
        decisionTreeReference} from '../controllers/db_calls.js'


const router : Router = express.Router()

router.get('/idcls', idCLs)
router.get('/tissues', tissues)
router.get('/experiments', experiments)
router.get('/databasestatistics', databaseStatistics)
router.get('/antibodiestable', antibodiesTable)

router.post('/findabs', findAbs)
router.post('/plotabs', plotAbs)
router.post('/findcelltypes', findCellTypes)
router.post('/plotcelltypes', plotCellTypes)
router.post('/findexperiments', findExperiments)
router.post('/whichantibodies', whichAntibodies)
router.post('/whichcelltypes', whichCellTypes)
router.post('/whichexperiments', whichExperiments)
router.post('/stveareference', stveaReference)
router.post('/convertcelltype', convertCellType)
router.post('/antibodypanelreference', antibodyPanelReference)
router.post('/decisiontreereference', decisionTreeReference)

// Website routes
router.post('/findabsweb', findAbsWeb)
router.post('/plotabsweb', plotAbsWeb)
router.post('/findcelltypesweb', findCellTypesWeb)
router.post('/plotcelltypesweb', plotCellTypesWeb)

// module.exports = router;
export default router