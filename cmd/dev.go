package main

import (
	"fmt"
	"math/rand"
	"time"

	"github.com/zedseven/bch"
)

func corruptData(data *[]uint8, errors int) {
	f := randIndex(len(*data))
	for i := 0; i < errors; i++ {
		l, _ := f()
		(*data)[l] ^= 1
	}
}

func randIndex(count int) func() (int, error) {
	pool := makeRange(int64(count))
	return func() (int, error) {
		if count <= 0 {
			return -1, nil
		}

		j := rand.Int63n(int64(count)) //I'm aware this isn't crypto/rand, but I needed to be able to seed it

		count--

		p := int(pool[j])

		pool[j] = pool[count]
		pool = pool[:count]

		return p, nil
	}
}

func makeRange(max int64) []int64 {
	r := make([]int64, max)
	for i := range r {
		r[i] = int64(i)
	}
	return r
}

func main() {
	rand.Seed(time.Now().Unix())


	/*for y := 1; y < 16; y++ {
		for x := 1; x < 16; x++ {
			res, _ := bch.TotalBitsForConfig(x, y)
			//res, _ := bch.TotalBitsForConfig(5, 8)
			fmt.Printf("%3d data bits, %2d error-fixes: %3d\n", x, y, res)
		}
	}*/


	data := []uint8{0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0}
	fmt.Printf("Original data: %v\n", data)
	config, err := bch.CreateConfig(265, 1)
	if err != nil {
		fmt.Println(err.Error())
	}
	code, err := bch.Encode(config, &data)
	if err != nil {
		fmt.Println(err.Error())
	}
	fmt.Printf("This is a %v.\n", config)
	fmt.Printf("Encoded data:   %v\n", code)
	//corruptData(&code, int(rand.Int63n(int64(12))))
	fmt.Printf("Corrupted data: %v\n", code)
	fmt.Println("Is data corrupt?", bch.IsDataCorrupted(config, code))
	recv, errors, err := bch.Decode(config, &code)
	if err != nil {
		fmt.Println(err.Error())
	}
	fmt.Printf("Decoded data:   %v\n", recv)
	fmt.Printf("There were %d error(s) in the corrupt data.", errors)

	fmt.Printf("\n\n")

	data = []uint8{0, 1, 0, 0, 1, 1, 0, 0}
	fmt.Printf("Original data: %v\n", data)
	code, err = bch.Encode(config, &data)
	if err != nil {
		fmt.Println(err.Error())
	}
	fmt.Printf("This is a %v.\n", config)
	fmt.Printf("Encoded data:   %v\n", code)
	corruptData(&code, int(rand.Int63n(int64(12))))
	fmt.Printf("Corrupted data: %v\n", code)
	fmt.Println("Is data corrupt?", bch.IsDataCorrupted(config, code))
	recv, errors, err = bch.Decode(config, &code)
	if err != nil {
		fmt.Println(err.Error())
	}
	fmt.Printf("Decoded data:   %v\n", recv)
	fmt.Printf("There were %d error(s) in the corrupt data.", errors)
}