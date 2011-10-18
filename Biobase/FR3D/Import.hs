{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}

-- | Importing of FR3D data. Both "basepairs" and "near interactions" are
-- currently supported. More parsers will come if required.

module Biobase.FR3D.Import where

import Control.Arrow
import Data.ByteString.Char8 as BS
import Data.Char
import Data.Iteratee as I
import Data.Iteratee.Char as I
import Data.Iteratee.IO as I
import Data.Iteratee.ListLike as I
import Data.List as L
import Data.Maybe
import qualified Data.Map as M
import System.FilePath.Find as F

import Biobase.Secondary

import Biobase.FR3D



-- | An Iteratee from a bytestring to one FR3D entry. Since each file contains
-- exactly one entry, this is no problem.

iFR3D :: (Monad m) => Iteratee ByteString m FR3D
iFR3D = joinI $ enumLinesBS f where
  f = do
    I.head -- fr3d header
    I.head -- sequence header
    cs' <- I.break ((/="#") . BS.take 1)
    I.head -- basepairs header
    xs <- stream2list -- and all basepairs
    let cs = L.map (second (BS.drop 1) . BS.span isAlphaNum . BS.drop 2) $ cs'
    return FR3D
      { pdbid = maybe "" (BS.take 4) $ listToMaybe xs
      , chains = cs
      , basepairs = {- L.map (fixSeqpos cs) . -} L.map bs2basepair $ xs
      }

{-
 - This would be for fixing sequence position information, but it seems that
 - FR3D does not store this info consistently...
 -
fixSeqpos :: [(ByteString,ByteString)] -> Basepair -> Basepair
fixSeqpos cs bp@Basepair{..} = bp{seqpos1 = seqpos1 - cl M.! chain1, seqpos2 = seqpos2 - cl M.! chain2} where
  cl = M.fromList . snd . L.mapAccumL f 0 $ cs
  f acc x = (acc + BS.length (snd x), (fst x, acc))
-}

-- | Helper function turning a bytestring line into a basepair entry

bs2basepair :: ByteString -> Basepair
bs2basepair s
  | L.length ws /= 10 = error $ "can't parse line: " ++ unpack s
  | otherwise = Basepair
    { interaction = threeChar . BS.unpack $ ws!!1
    , nucleotide1 = BS.head $ ws!!2
    , pdbnumber1  = maybe (-1) fst . readInt $ ws!!3
    , chain1      = ws!!4
    , seqpos1     = maybe (-1) (subtract 1 . fst) . readInt $ ws!!5
    , nucleotide2 = BS.head $ ws!!6
    , pdbnumber2  = maybe (-1) fst . readInt $ ws!!7
    , chain2      = ws!!8
    , seqpos2     = maybe (-1) (subtract 1 . fst) . readInt $ ws!!9
    }
  where ws = BS.words s

-- | Convenience function: given a directory name, extracts a list of all FR3D
-- entries.

fromDirSelect :: String -> FilePath -> IO [FR3D]
fromDirSelect select fp = do
  fs <- F.find always (fileName ~~? select) fp
  mapM (\f -> run =<< enumFile 8192 f iFR3D) fs

-- | This one select the "near interactions"

fromDirNear = fromDirSelect "*near_interactions_FR3D.txt"

-- | And this one the "basepairs" (this one you normally want).

fromDir = fromDirSelect "*basepairs_FR3D.txt"
