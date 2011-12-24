{-# LANGUAGE BangPatterns #-}
{-# OPTIONS_GHC -O2 -fexcess-precision -funbox-strict-fields #-}

import qualified Data.ByteString.Char8 as B
import Data.Function (on)
import Data.List (foldl', sortBy, transpose)
import qualified Data.IntSet as S
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import System.Environment
import Text.Printf


type ValueClassMarginal = (U.Vector Int, [(Int, Double)])


mutualInfoInnerLoop :: Double -> U.Vector (Int, Int) -> Double -> (Int, Int, Double) -> Double
mutualInfoInnerLoop n xys !acc (!i, !j, !px_py)
    | px_py == 0 || pxy == 0 = acc
    | otherwise              = pxy * logBase 2 ( pxy / px_py ) + acc
    where
        pxy = ({-# SCC "foldr'" #-} fromIntegral . U.foldr' (accumEq2 i j) 0 $ xys ) / n
        accumEq2 :: Int -> Int -> (Int, Int) -> Int -> Int
        accumEq2 !i !j (!i', !j') !acc
            | i' == i && j' == j = acc + 1
            | otherwise          = acc


mutualInfo :: Double -> ValueClassMarginal -> ValueClassMarginal -> Double
mutualInfo n (xs, xcp) (ys, ycp) = foldl' ({-# SCC "mutualInfoInnerLoop" #-} mutualInfoInnerLoop n $ U.zip xs ys) 0 [ (i, j, px * py) | (!i, !px) <- xcp, (!j, !py) <- ycp ]


maxRel :: Double -> V.Vector ValueClassMarginal -> ValueClassMarginal -> U.Vector Double
maxRel n xcls ycls = U.generate (V.length xcls) (\i -> mutualInfo n (xcls V.! i) ycls)


data MrmrMethod = MID | MIQ deriving (Eq, Show)


doMrmr :: Int -> MrmrMethod -> V.Vector ( U.Vector Int ) -> U.Vector Int -> (U.Vector (Int, Double), U.Vector (Int, Double))
doMrmr k method cols ys
    | n == 0    = (U.empty, U.empty)
    | otherwise = (imaxrels, mrmrRecurser (k - 1) method n xcls imaxrels $ U.singleton . U.maximumBy (compare `on` snd) $ imaxrels)
    where
        n               = fromIntegral . U.length $ ys
        imarginals :: U.Vector Int -> [(Int, Double)]
        imarginals vals = map (\i -> (i, fromIntegral ( U.foldl' (\l r -> if r == i then l + 1 else l) 0 vals ) / n)) $ S.elems . S.fromList . U.toList $ vals
        ycls            = (ys, imarginals ys)
        xcls            = V.zip cols $ V.map imarginals cols
        imaxrels        = U.zip (U.enumFromN 0 $ U.length maxrels) maxrels
            where
                maxrels = maxRel n xcls ycls


mrmrRecurser :: Int -> MrmrMethod -> Double -> V.Vector ValueClassMarginal -> U.Vector (Int, Double) -> U.Vector (Int, Double) -> U.Vector (Int, Double)
mrmrRecurser k method n xcls maxrels mrmrs
    | k == 0    = mrmrs
    | otherwise = mrmrRecurser (k - 1) method n xcls maxrels $ U.snoc mrmrs $ mrmrNextIdx method n xcls maxrels mrmrs


mrmrNextIdx :: MrmrMethod -> Double -> V.Vector ValueClassMarginal -> U.Vector (Int, Double) -> U.Vector (Int, Double) -> (Int, Double)
mrmrNextIdx method n xcls maxrels mrmrs =
    U.foldl' folder (0, 0) maxrels
    where
        folder l@(_, v) (i, maxrel)
            | U.elem i ( U.map fst mrmrs ) || v' <= v = l
            | otherwise                               = (i, v')
            where
                v' = mrmrVal method n xcls i maxrel mrmrs


mrmrVal :: MrmrMethod -> Double -> V.Vector ValueClassMarginal -> Int -> Double -> U.Vector (Int, Double) -> Double
mrmrVal method n xcls i maxrel mrmrs
    | method == MID = maxrel - mrmr
    | method == MIQ = maxrel / mrmr
    where
        mrmr = U.foldl' (\l r -> mutualInfo n ( xcls V.! fst r ) ( xcls V.! i ) + l) 0 mrmrs / fromIntegral ( U.length mrmrs )


parseCsv :: B.ByteString -> ([String], V.Vector ( U.Vector Int ), U.Vector Int)
parseCsv input = (map B.unpack $ B.split ',' labels, cols, ys)
    where
        cols           = V.fromList coldata
        ys:coldata     = map U.fromList . transpose $ [ map readInt $ B.split ',' row | row <- rowdata ]
        labels:rowdata = [ line | line <- B.lines input, line /= B.empty ]


readInt :: B.ByteString -> Int
readInt str = case B.readInt str of
    Just (v, _) -> v
    Nothing     -> error "non-integer value is unparseable"


main = do
        [ks, f] <- getArgs
        let k = read ks :: Int
        s <- B.readFile f
        let (labels, cols, ys) = parseCsv s
        let (maxrels, mrmrs) = doMrmr k MID cols ys
        printf "*** MaxRel features ***\n"
        sequence_ [ printf "%d\t%s\t%.3f\n" ( idx + 1 ) ( labels !! ( idx + 1 ) ) maxrel | (idx, maxrel) <- take k . reverse . sortBy (compare `on` snd) . U.toList $ maxrels ]
        printf "\n*** mRMR features ***\n"
        sequence_ [ printf "%d\t%s\t%.3f\n" ( idx + 1 ) ( labels !! ( idx + 1 ) ) mrmr   | (idx, mrmr)   <- U.toList mrmrs   ]
