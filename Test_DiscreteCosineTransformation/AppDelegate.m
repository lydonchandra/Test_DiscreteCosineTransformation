//
//  AppDelegate.m
//  Test_DiscreteCosineTransformation
//
//  Created by Don on 10/21/15.
//  Copyright © 2015 Don. All rights reserved.
//

#import "AppDelegate.h"
#include "dct.h"
#include "targa.h"
#include "dct_2d.h"

@interface AppDelegate ()

@end

@implementation AppDelegate

extern uint64_t dispatch_benchmark(size_t count, void (^block)(void));


void perf_test_dct_8_elements() {
    static size_t const iterations = 200000;
    
    const int COUNT = 8;
    
    double *out_1D_array = malloc(COUNT* sizeof(double));
    
    //count = 8
    //8.0, 16.0, 24.0, 32., 40., 48., 56., 64.
    //dct_ii: 101.823376 -51.538584 -0.000000 -5.387638 0.000000 -1.607223 -0.000000 -0.405619
    //llm_dct: 101.823376 -51.538584 0.000000 -5.387638 0.000000 -1.607223 0.000000 -0.405619
//    2015-10-22 23:32:47.444 Test_DiscreteCosineTransformation[741:437739] llm_dct Avg. Runtime: 401 nanoseconds
//    2015-10-22 23:32:47.551 Test_DiscreteCosineTransformation[741:437739] aan_dct Avg. Runtime: 525 nanoseconds
//    2015-10-22 23:32:48.404 Test_DiscreteCosineTransformation[741:437739] dct_ii Avg. Runtime: 4262 nanoseconds
    
    //idct_ii: 8.486029 15.099155 24.325407 32.469266 39.530734 47.674593 56.900845 63.513971
//    NSLog(@"%f %f %f %f %f %f %f %f",
//          out_1D_array[0],out_1D_array[1], out_1D_array[2], out_1D_array[3],
//          out_1D_array[4], out_1D_array[5], out_1D_array[6], out_1D_array[7]);
    
    uint64_t average_runtime = dispatch_benchmark(iterations, ^{
        double in_1D_array[COUNT] = {
            8., 16., 24., 32., 40., 48., 56., 64.
        };
        llm_dct(in_1D_array, out_1D_array);
    });
    NSLog(@"llm_dct Avg. Runtime: %llu nanoseconds", average_runtime);

    average_runtime = dispatch_benchmark(iterations, ^{
        double in_1D_array[COUNT] = {
            8., 16., 24., 32., 40., 48., 56., 64.
        };
        aan_dct(in_1D_array, out_1D_array);
    });
    NSLog(@"aan_dct Avg. Runtime: %llu nanoseconds", average_runtime);

    average_runtime = dispatch_benchmark(iterations, ^{
        double in_1D_array[COUNT] = {
            8., 16., 24., 32., 40., 48., 56., 64.
        };
        dct_ii_don(COUNT, in_1D_array, out_1D_array);
    });
    NSLog(@"dct_ii Avg. Runtime: %llu nanoseconds", average_runtime);
    
}


void perf_test_targa() {
    NSString *inTgaFilePath = [[NSBundle mainBundle] pathForResource:@"in" ofType:@"tga"];
    
    
    NSArray  *dirPaths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
    NSString *docsDir = [dirPaths objectAtIndex:0];
    
    NSDateFormatter *formatter = [[NSDateFormatter alloc] init];
    [formatter setDateFormat:@"yyyy_MM_dd_HH:mm:ss.SSS"];
    
    NSDate   *currentDate = [NSDate date];
    
    NSString *currentDateString = [formatter stringFromDate:currentDate];
    
    NSString *filename = [NSString stringWithFormat:@"raw%@.tga", currentDateString];
    NSString *databasePath = [[NSString alloc] initWithString:[docsDir stringByAppendingPathComponent:filename]];
    NSString *outTgaFilePath = databasePath;
    const char *outTgaFilePath_c = [outTgaFilePath UTF8String];
    char * outTgaFilePath_c_copy = calloc([databasePath length]+1, 1);
    strncpy(outTgaFilePath_c_copy, outTgaFilePath_c, [outTgaFilePath length]);
    
    const char *inTgaFilePath_c = [inTgaFilePath UTF8String];
    char *cpy = calloc([inTgaFilePath length]+1, 1);
    strncpy(cpy, inTgaFilePath_c, [inTgaFilePath length]);
    
    const size_t iterations = 30;
    uint64_t average_runtime = dispatch_benchmark(iterations, ^{
        main_old2(cpy, outTgaFilePath_c_copy);        
    });
    NSLog(@"perf test targa dct=>quantize=>idct: %llu milliseconds", average_runtime/1000000);
//    main_old2(cpy, outTgaFilePath_c_copy);
    
    UIImage *savedImage =  [[UIImage alloc] initWithContentsOfFile:databasePath];
    UIImageWriteToSavedPhotosAlbum(savedImage, nil, nil, nil);

}



void perf_test_dct() {
    static size_t const iterations = 200000;
    
    const int COUNT = 16;
    
    double *out_1D_array = malloc(COUNT* sizeof(double));
    
    uint64_t average_runtime = dispatch_benchmark(iterations, ^{
        double in_1D_array[COUNT] = {
            8., 16., 24., 32., 40., 48., 56., 64., 72., 80.
            , 88., 96., 104., 112., 120., 128.
        };
        dct_ii_don(COUNT, in_1D_array, out_1D_array);
    });
    NSLog(@"dct Avg. Runtime: %llu nanoseconds", average_runtime);
    
    average_runtime = dispatch_benchmark(iterations, ^{
        double in_1D_array[COUNT] = {
            272.000000,-146.492248,-0.000000,-16.060225,-0.000000,-5.612698
            ,-0.000000,-2.716334,-0.000000,-1.501422,-0.000000,-0.857121
            ,-0.000000, -0.448301,0.000000,-0.139962
        };
        idct_ii_don(COUNT, in_1D_array, out_1D_array);
    });
    NSLog(@"inverse dct Avg. Runtime: %llu nanoseconds", average_runtime);

    
    
    for (int idx=0; idx<COUNT; idx++) {
//        printf("%f ", out_1D_array[idx]);
    }
//    printf("\n");
    
    //zeroing idx 5,6,7
    out_1D_array[5] = out_1D_array[6] = out_1D_array[7] = 0.;
    
    double *out_1D_array_inverse = malloc(COUNT * sizeof(double));
    idct_ii_don(COUNT, out_1D_array, out_1D_array_inverse);
    for (int idx=0; idx<COUNT; idx++) {
//        printf("%f ", out_1D_array_inverse[idx]);
    }
//    printf("\n");    
}

- (BOOL)application:(UIApplication *)application didFinishLaunchingWithOptions:(NSDictionary *)launchOptions {
    
    perf_test_targa();
    //perf_test_dct_8_elements();
    //perf_test_dct();
    
    
    // Override point for customization after application launch.
    return YES;
}

- (void)applicationWillResignActive:(UIApplication *)application {
    // Sent when the application is about to move from active to inactive state. This can occur for certain types of temporary interruptions (such as an incoming phone call or SMS message) or when the user quits the application and it begins the transition to the background state.
    // Use this method to pause ongoing tasks, disable timers, and throttle down OpenGL ES frame rates. Games should use this method to pause the game.
}

- (void)applicationDidEnterBackground:(UIApplication *)application {
    // Use this method to release shared resources, save user data, invalidate timers, and store enough application state information to restore your application to its current state in case it is terminated later.
    // If your application supports background execution, this method is called instead of applicationWillTerminate: when the user quits.
}

- (void)applicationWillEnterForeground:(UIApplication *)application {
    // Called as part of the transition from the background to the inactive state; here you can undo many of the changes made on entering the background.
}

- (void)applicationDidBecomeActive:(UIApplication *)application {
    // Restart any tasks that were paused (or not yet started) while the application was inactive. If the application was previously in the background, optionally refresh the user interface.
}

- (void)applicationWillTerminate:(UIApplication *)application {
    // Called when the application is about to terminate. Save data if appropriate. See also applicationDidEnterBackground:.
}

@end
