import express, { Express, Router, Request, Response } from 'express'
import { idCLs, runLMM } from '../controllers/db_calls.js'


const router : Router = express.Router()

router.get('/idCLs', idCLs)
router.post('/abLMM', runLMM)

// module.exports = router;
export default router