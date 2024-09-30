#include "cut_rf_up.h"
#include "liver_cut.h"  // 包含你已经模块化的函数，比如 plane_generator, cut_run2 等
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第六刀的操作
vtkSmartPointer<vtkPolyData> cut_rf_up(vtkSmartPointer<vtkPolyData> liver_r2, 
               vtkSmartPointer<vtkPolyData> rb1v, 
               vtkSmartPointer<vtkPolyData> rf_down_v, 
               vtkSmartPointer<vtkPolyData> rf_up_v,
              double (&mbounds)[6]) 
{
    // 第6刀，右前上段
    auto exc_rf_up_v = mergePolyData(rf_down_v, rb1v);
    auto [rf_up_b, rf_up_v_points] = find_points(exc_rf_up_v, rf_up_v);
    auto imp6 = plane_generator(rf_up_b);
    auto [liver_r3, rf_up_v_poly] = cut_run2(imp6, liver_r2,mbounds);

    // 保存右前上段切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/rf_up.stl");
    stlWriter->SetInputData(rf_up_v_poly);
    stlWriter->Write();

    vtkSmartPointer<vtkSampleFunction> sampleFunction = vtkSmartPointer<vtkSampleFunction>::New();
    sampleFunction->SetImplicitFunction(imp6); // 您的隐式函数
    sampleFunction->SetModelBounds(mbounds[0], mbounds[1], mbounds[2], mbounds[3], mbounds[4], mbounds[5]); // 设置采样范围
    sampleFunction->SetSampleDimensions(50, 50, 50); // 设置采样分辨率
    sampleFunction->ComputeNormalsOff();
    sampleFunction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(sampleFunction->GetOutputPort());
    contourFilter->SetValue(0, 0.0); // 提取等值为 0 的面，即隐式函数的零等值面
    contourFilter->Update();

    stlWriter->SetFileName("C:/code/liver_cut/new_imp/imp_rf_up.stl"); // 设置输出文件名
    stlWriter->SetInputConnection(contourFilter->GetOutputPort());
    stlWriter->Write();


    std::cout << "Successfully saved rf_up cut" << std::endl;
    std::cout << "rf_up cut completed" << std::endl;

    return liver_r3;
}
