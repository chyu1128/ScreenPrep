import pytest
from unittest.mock import MagicMock
from app import ScreenPrep, Record

def test_record_dataclass():
    """驗證 Record 資料結構是否正常初始化"""
    rec = Record(compound_id="CMP1", original_smiles="CCO", value="2.5")
    assert rec.compound_id == "CMP1"
    assert rec.outcome == "Pending"

def test_screen_prep_initial_state(monkeypatch):
    """測試 App 初始化時的初始狀態（使用 monkeypatch 避免 Tkinter 真正跳出視窗）"""
    # 阻斷 mainloop 避免測試卡住
    monkeypatch.setattr(ScreenPrep, "mainloop", lambda self: None)
    
    app = ScreenPrep()
    assert len(app.records) == 0
    assert app.processed is False
    # 關閉視窗以釋放資源
    app.destroy()

def test_process_logic_with_mock_data(monkeypatch):
    """測試核心的預處理邏輯（包含合法、重複與不合法的 SMILES）"""
    monkeypatch.setattr(ScreenPrep, "mainloop", lambda self: None)
    
    app = ScreenPrep()
    # 模擬注入測試資料
    app.records = [
        Record("ID1", "CCO", "1.2"),        # 正常乙醇分子
        Record("ID2", "CCO.Cl", "3.4"),     # 帶有鹽酸鹽的乙醇分子（脫鹽後應與第一個重複）
        Record("ID3", "INVALID_SMILES", "5.6") # 不合法的結構
    ]
    
    # 模擬 RDKit 處理時彈出的提示視窗，避免測試中斷
    app._show_running_dialog = MagicMock()
    
    # 執行預處理
    app.process()
    
    # 驗證輸出結果
    # 1. 第一個正常分子應該進入 Ready 狀態
    assert app.records[0].outcome == "Ready"
    assert app.records[0].cleaned_smiles == "CCO"
    
    # 2. 第二個分子去鹽後變成 CCO，與第一個重複，應被歸類為 Excluded 重複項
    assert app.records[1].outcome == "Excluded"
    assert "Duplicate" in app.records[1].detail
    
    # 3. 第三個分子是不合法 SMILES，應被排除
    assert app.records[2].outcome == "Excluded"
    assert "Invalid SMILES" in app.records[2].detail
    
    app.destroy()