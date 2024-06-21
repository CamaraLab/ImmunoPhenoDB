import express, { Express, Router, Request, Response } from 'express'
import {idCLs, 
        tissues,
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
        antibodyPanelReference } from '../controllers/db_calls.js'


const router : Router = express.Router()

router.get('/idcls', idCLs)
router.get('/tissues', tissues)
router.get('/databasestatistics', databaseStatistics)

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

// module.exports = router;
export default router